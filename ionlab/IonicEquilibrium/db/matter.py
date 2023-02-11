"""
Provide classes for working with all kinds of chemistry substances, mixtures
"""
from ionlab.IonicEquilibrium.db.dbmanage import loc
import sympy
from ionlab.conf import IonicEquilibriumOptions
from typing import Union
import copy


class InitIonQuantity:
    def __init__(self, M: float = 0, mol: float = 0, mass: float = 0):
        self.M = M
        self.mol = mol
        self.mass = mass

    def add(self, values):
        self.M = self.M + values.get('M', 0)
        self.mol = self.mol + values.get('mol', 0)
        self.mass = self.mass + values.get('mass', 0)

    def update(self, values):  # Like a dictionary update
        self.M = values.get('M', self.M)
        self.mol = values.get('mol', self.mol)
        self.mass = values.get('mass', self.mass)

    def recalculate(self, volume, mol: float = None):
        if mol is not None:
            self.mol = mol
        self.M = self.mol/volume

    def __add__(self, other):
        self.M = self.M + other.M
        self.mol = self.mol + other.mol
        self.mass = self.mass + other.mass
        return self

    def __getitem__(self, item):
        return self.__dict__[item]

    def __repr__(self):
        return str(self.__dict__)

    def __str__(self):
        return str(self.__dict__)


class Ion:
    def __init__(self, formula: str, chemtype: str, count_in_compound: int = None, datadict=None):
        if chemtype not in ('cation', 'anion'):
            raise TypeError('Invalid chemtype for ion initialization.')
        self.count = count_in_compound
        if datadict is None:
            ion_datadict = loc(formula)
        else:
            ion_datadict = datadict
        self.datadict = ion_datadict
        self.formula = ion_datadict['formula']
        self.name = ion_datadict['name']
        self.pKa = ion_datadict['pKa']
        self.species = ion_datadict['species']
        self.conjugate_count = ion_datadict['conjugate_count']
        self.charge = ion_datadict['charge']
        self.multiplier = ion_datadict['multiplier']
        self.init = InitIonQuantity(0, 0)

    def __getitem__(self, item):
        return self.__dict__[item]


class SpecialIon(Ion):
    """You shouldn't initialize this class yourself. This is called from matter.Anion when it detects a special ion"""
    def __init__(self, formula: str, tier: int, pKspec: float, relation: str):
        if tier < 1:
            msg = 'pKspec numbering must start from 1'
            raise ValueError(msg)
        super().__init__(formula, 'anion')
        self.formula = formula
        self.pKspec = {tier: pKspec}
        self.relation = relation
        self.tier = tier

    def __repr__(self):
        special_details = f'matter.SpecialIon: {self.formula}, pKspec[{self.tier}] = {self.pKspec[self.tier]},\
{self.relation}'
        return special_details

    def __getitem__(self, item):
        return self.__dict__[item]


class Cation(Ion):
    def __init__(self, formula: str = None, count_in_compound: int = None, datadict=None):
        super().__init__(formula, 'cation', count_in_compound, datadict)

    def __repr__(self):
        cation_details = f'matter.Cation: {self.formula}, count = {self.count}, {self.init.M}M, {self.init.mol} mol'
        return cation_details


class Anion(Ion):
    def __init__(self, formula: str = None, count_in_compound: int = None, datadict=None):
        super().__init__(formula, 'anion', count_in_compound, datadict)
        self.proton_count = self.datadict.get('proton_count', 0)
        self.specials = self.datadict.get('special', None)
        if self.specials is not None:
            temp = self.specials
            self.specials = []
            for i, spec in enumerate(temp):
                special_ion = SpecialIon(spec['formula'], i + 1, spec['pKspec'][i + 1], spec['relation'])
                self.specials.append(special_ion)

    def __repr__(self):
        anion_details = f'matter.Anion: {self.formula}, count = {self.count}, {self.init.M}M, {self.init.mol} mol'
        return anion_details


class Compound:
    def __init__(self, *args, **kwargs):
        """
        Smartly construct a compound from 2 methods: Inputting the formula directly or the cation, anion and their count
        Also supports dictionary input for ions, as long as it's the only input
        :param args: [formula] or in order [cation, anion], [cation_count, anion_count]
            :note The order of cations and anions doesn't matter. However, cation_count and anion_count should be in
            the same order as the cations and anions and you MUST NOT mix args and kwargs together. Use one or another
                Compound('Ba2+', 'SO42-', 1, 1)
                Compound(1 (SO42-), 1 (Ba2+), 'SO42-', 'Ba2+')
                Compound('BaSO4')
                Compound({'Ba2+': 1, 'SO42-': 1})
                Compound(['Ba2+', 'SO42-'], [1, 1]) TODO: This
        :param kwargs: [formula] or [cation, anion, cation_count, anion_count]
        """

        self.options = IonicEquilibriumOptions
        self.net_charge = -1
        formula, cations, anions, cation_count, anion_count = self._validate_creation(args, kwargs)  # Validate input

        self.cations = list()
        self.anions = list()

        # METHOD 1: Given the formula
        if formula is not None:
            datadict = loc(formula)
            # Special cases: Bases that are electrically neutral. They shall be treated as an anion
            self.type = datadict['type']
            self.formula = datadict['formula']
            if self.type == 'acid':  # For acid
                self.init_cation(cation='H+', cation_count=datadict['conjugate_count'])
                conjugate_base_formula = datadict['species'][-1]
                self.init_anion(anion=conjugate_base_formula, anion_count=1)
            if self.type == 'base':  # For base like NH3, NH2OH
                self.init_cation(cation='H+', cation_count=0)
                self.init_anion(anion=formula, anion_count=1)

        else:  # METHOD 2: Cations and anions
            for i, cation in enumerate(cations):
                self.init_cation(cation, cation_count[i])

            # Check if there are multiple anion and add their info
            for i, anion in enumerate(anions):
                self.init_anion(anion, anion_count[i])

        # Checks that the compound created satisfy charge conservation rule
        self.ensure_neutral_charge()

    @staticmethod
    def _validate_creation(args: tuple, kwargs: dict) -> (str, str, str, int, int):
        """
        Validate, type check arguments and extract the creation method and details.
        :param args: args from __init__
        :param kwargs: kwargs from __init__
        :return: formula, cation, anion, cation_count, anion_count
        """
        formula, cation, anion, cation_count, anion_count = \
            kwargs.get('formula'), kwargs.get('cation', []), kwargs.get('anion', []), \
            kwargs.get('cation_count', []), kwargs.get('anion_count', [])
        if not isinstance(cation, list):
            cation = [cation]
        if not isinstance(anion, list):
            anion = [anion]
        if not isinstance(cation_count, list):
            cation_count = [cation_count]
        if not isinstance(anion_count, list):
            anion_count = [anion_count]
        ions = []
        counts = []

        for i, arg in enumerate(args):  # Extract formula, and cation, anion, into ions and counts for unpacking
            if isinstance(arg, str):
                if '+' in arg or '-' in arg:
                    ions.append(arg)
                else:
                    formula = arg
            if isinstance(arg, int):
                counts.append(arg)
        for ion, count in zip(ions, counts):  # Extract cation_count, anion_count using the cation, anion ordering
            if '+' in ion:
                cation.append(ion)
                cation_count.append(count)
            if '-' in ion:
                anion.append(ion)
                anion_count.append(count)

        # Dictionary input
        if len(args) == 1 and isinstance(args[0], dict):
            ions = list(args[0].keys())
            cation.append(ions[0] if '+' in ions[0] else ions[1])
            anion.append(ions[0] if '-' in ions[0] else ions[1])
            cation_count.append(args[0][cation])
            anion_count.append(args[0][anion])

        msg1 = f'{formula=}, {cation=}, {anion=}, {cation_count=}, {anion_count=}'
        if formula is None and not all([cation, anion, cation_count, anion_count]):  # No formula, no ions
            msg2 = 'Not sufficient data to create a compound.'
            raise ValueError(msg2 + '\n' + msg1)
        if formula is not None and all([cation, anion, cation_count, anion_count]):  # Formula and ion
            msg2 = 'Too much data was give. Either input the formula or the ions'
            raise ValueError(msg2 + '\n' + msg1)
        return formula, cation, anion, cation_count, anion_count

    def init_cation(self, cation, cation_count):
        cation = Cation(cation, cation_count)
        self.cations.append(cation)

    def init_anion(self, anion, anion_count):
        anion = Anion(anion, anion_count)
        if anion.proton_count > 0:  # If there are still free protons in the anion, break it down into H+ and base.
            self.init_cation('H+', cation_count=anion.proton_count)
            anion = Anion(anion.species[-1], anion_count)
        self.anions.append(anion)

    def ensure_neutral_charge(self):
        net_charge = 0
        for cation in self.cations:
            net_charge += cation.charge * cation.count
        for anion in self.anions:
            net_charge += anion.charge * anion.count
        e = f'The compound created is not a neutral compound as net charge = {net_charge}'
        self.net_charge = net_charge
        if self.net_charge != 0:
            raise AttributeError(e)

    def __getitem__(self, item):
        return self.__dict__[item]


class Solution:
    def __init__(self, volume):
        """
        :param volume: In terms of litres (cubic decimeter)
        """
        if not isinstance(volume, int) and not isinstance(volume, float):
            raise TypeError('Invalid type for volume.')
        self._volume = volume
        self.anions = {}
        self.cations = {}
        self._pH = None
        self.options = IonicEquilibriumOptions

    @property
    def volume(self):
        return self._volume

    @volume.setter
    def volume(self, volume: Union[float, int]):
        self._volume = volume
        for formula, anion in self.anions.items():
            anion.init.recalculate(volume)
        for formula, cation in self.cations.items():
            cation.init.recalculate(volume)

    def add_ion(self, ion: Union[Cation, Anion], init_mol: float = 0) -> None:
        if not isinstance(ion, Cation) and not isinstance(ion, Anion):
            raise TypeError(f'Invalid type for ion: {type(ion)}')
        """Adding ions only requires its init_mol. Init_M will automatically be calculated at the last stage."""
        init_M = sympy.N(init_mol / self.volume, self.options.PRECISION)
        new_init_quantity = InitIonQuantity(init_M, init_mol)
        # New ions will be initialized by the Solution object, so it will not inherit init quantity.
        ion = copy.deepcopy(ion)
        ion.init.update({'mol': 0, 'M': 0})  # Unnecessary, but make sure that the Ion object's quantity is resetted.
        if isinstance(ion, Cation):  # For adding cation
            self.cations[ion.formula] = ion
            self.cations[ion.formula].init += new_init_quantity

        if isinstance(ion, Anion):  # For adding anion
            self.anions[ion.formula] = ion
            self.anions[ion.formula].init += new_init_quantity

    def add(self, sol: Compound, init_M: float = 0, init_mol: float = 0) -> None:
        """
        Add a solute (compound) to the solution (water is the default solvent)
        :param sol: A Compound object
        :param init_M: Molarity of the compound (M)
        :param init_mol: Mole of the compound (mol)
        :return:
        """
        if init_mol == 0 and init_M == 0:
            raise ValueError('init_M and init_mol can\'t both be 0.')
        init_mol = sympy.N(init_M * self.volume, self.options.PRECISION) if init_mol == 0 else init_mol
        for anion in sol.anions:
            self_mol = self.anions.get(anion.formula).init.mol if self.anions.get(anion.formula) is not None else 0
            self.add_ion(anion, init_mol=init_mol * anion.count + self_mol)
        for cation in sol.cations:
            self_mol = self.cations.get(cation.formula).init.mol if self.cations.get(cation.formula) is not None else 0
            self.add_ion(cation, init_mol=init_mol * cation.count + self_mol)

    def __add__(self, other):
        new_volume = self.volume + other.volume
        self.volume = new_volume
        # Update ions that exist in other
        for formula, anion in other.anions.items():
            other_mol = anion.init.mol
            # Reset other's ion quantity (so it can be added in case it's new)
            anion.init = InitIonQuantity()
            self_mol = self.anions[formula].init.mol if self.anions.get(formula) is not None else 0
            self.add_ion(anion, sympy.N(self_mol + other_mol, self.options.PRECISION))  # So we can calculate from the beginning

        for formula, cation in other.cations.items():
            other_mol = cation.init.mol
            # Reset other's ion quantity (so it can be added in case it's new)
            cation.init = InitIonQuantity()
            self_mol = self.cations[formula].init.mol if self.cations.get(formula) is not None else 0
            self.add_ion(cation, sympy.N(self_mol + other_mol, self.options.PRECISION))
        return self

    def __repr__(self):
        return str(self.__dict__)
