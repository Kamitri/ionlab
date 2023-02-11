from __future__ import annotations
from ionlab.IonicEquilibrium.db import matter
from ionlab.IonicEquilibrium.db.dbmanage import loc, iloc
import sympy
from sympy.core.relational import Equality
from ionlab.IonicEquilibrium.solver.gekko_solver import GEKKOSolver
from ionlab.IonicEquilibrium.solver.wolfram_solver import WolframSolver
from ionlab.conf import IonicEquilibriumOptions
from ionlab.IonicEquilibrium.utils import extension
from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wl, wlexpr
from itertools import product
import matplotlib.pyplot as plt
import numpy as np
import warnings


# TODO: Ion molarity = 0, delete ion
class Lab:
    def __init__(self, solution: matter.Solution):
        self.SOLVERS = {'wolfram': self.pH_solve_wolfram, 'gekko': self.pH_solve_gekko}
        self.solution = solution
        # Quick reference to all solutes
        self.anions = solution.anions
        self.cations = solution.cations

        # Reference dictionary used for keeping track of analysis
        self.constants = dict()
        self.init_anion_concentration = dict()
        self.init_cation_concentration = dict()
        self.anion_reference = dict()
        self.cation_reference = dict()
        self.complex_reference = dict()
        self.solid_reference = dict()
        self.logbook = {'H+': 'h'}
        self.calculated_eq = False
        # All equations used for referencing and debugging
        self.equations = dict()
        # Options
        self.options = IonicEquilibriumOptions

    def titration(self, titrant: matter.Compound, init_M: float | int,  max_vol: float | int,
                  step: float | int = 10**-3, solver: str = 'wolfram', plt_backend: str = 'qtagg') -> plt.figure:
        pHs = list()
        for curr_vol in np.arange(10**-5, max_vol, step):
            sol2 = matter.Solution(curr_vol)
            sol2.add(titrant, init_M=0.1)
            self.solution = self.solution + sol2
            lab = Lab(self.solution)
            lab.analyze_eq(solver=solver)
            pH = lab.get_pH()
            pHs.append(pH)
        plt.switch_backend(plt_backend)
        plt.plot(list(np.arange(10**-5, max_vol, step)), pHs)
        plt.xlabel(f'Titration Volume (mL)')
        plt.ylabel('pH')
        plt.title('Titration curve')
        return plt

    def get_pH(self, solver: str = 'gekko'):
        if not self.calculated_eq:
            self.analyze_eq(solver=solver)
        return self.eq_molarity_dict['pH']

    def analyze_eq(self, solver: str = 'gekko'):
        """
        One more fight. One more night.
        Analysis method:
            We will create 2 types of equation
                1) Relation between the most characteristic ion of a species with its init_C, constants and other ion members
                2) Ion conservation equation
        :param (str) solver: Available solvers are 'gekko'
        :return:
        """
        # To begin analysis, we create essentials dictionaries with useful infos to perform analysis
        self._preprocess()
        # Next we create equations relating the constant and other ions and update self.constants
        self._construct_const_eq()
        # Next we create equations relating initial concentration and other ions and
        # update self.init_concentration (both)
        self._construct_init_M_eq()
        # Next we create equations that can calculate the concentration of any ions at equilibrium and update
        self._construct_ion_value_eq()
        # Next we create solid equations (if any)
        self._construct_solid_eq()
        # Create ion conservation equations
        self._construct_ion_conservation_eq()

        # Finally, solve the equations
        self.eq_molarity_dict = self.SOLVERS[solver]()
        self._update_reference_at_eq()
        self.calculated_eq = True
        return self.eq_molarity_dict

    def _preprocess(self):
        """Construct 7 dictionary of references.
            1) anion_reference is for storing symbol of each anion in equilibrium, charge, formula, its equation and
                relationship...
            2) cation_reference is for storing symbol of each cation in equilibrium, charge, formula, its equation and
                relationship...
            3) complex_reference is for storing symbol of each complex in equilibrium, charge, formula, its equation and
                relationship...
            4) solid_reference is for storing symbol of each precipitation in equilibrium, charge, formula, its equation and
                relationship...
            5) init_anion_concentration is for storing initial concentration of anion species (e.g: PO₄³⁻, HPO₄²⁻)
            6) init_cation_concentration is for storing initial concentration of cation species (e.g: Na⁺, Ba²⁺)
            7) constants is for storing constants like pKa, pKw or pβ
        1) (dict) anion_reference / (dict) cation_reference
            [key] structure:
                'x' + 'ani'(anion) / 'cat'(cation) + {id} + 'st' + {state}(conjugated_state)
                'h' (for proton H⁺)
            [value] structure:
                'formula': (str) The formula of that ion (e.g: 'PO₄³⁻')
                'id': (int) The cation / anion id of that ion species (start counting from 0, cation and anion id are different)
                'charge': (int) Net charge of that ion (positive for cation and negative for anion)
                'state': (int) Conjugated state of that ion (0 is conjugated acid, n is conjugated base)
                'eq': (sympy.Eq) An equation relating this ion with constants and other ions
                'value': (sympy.expr | float) Can either be an expression of C at equilibrium of this ion or C itself
                'multiplier': (int) How many element of this species is contained within one ion (usually 1; 2 for HF2-...)
                'const_symbol': (str) The constant related to this ion (state 0 has no constant)
                'origin': (str) How this ion is related to other ion ['acid-base', 'special']

        2) (dict) complex_reference
            [key] structure:
                'x' + 'cat' + {cation_id} + 'ani' + {anion_id} + 'cp' + {state}(complex_state)
            [value] structure:
                'formula': (str) The formula of this complex (e.g: 'AgBr2-')
                'cation': (str) The formula of the cation of this complex
                'ligand': (str) The formula of the anion of this complex
                'cation_symbol': (str) The symbol of the cation of this complex, for easier access
                'ligand_symbol': (str) The symbol of the ligand of this complex, for easier access
                'cat_id': (int) The id of the cation of this complex
                'ani_id': (int)The id of the anion of this complex
                'charge': (int) Net charge of this complex
                'cation_count': (int) How many cation is in one complex ion
                'ligand_count': (int) How many ligand is in one complex ion
                'state': (int) Complex state of this complex (1 is the lowest, n is the highest)
                'eq': (sympy.Eq) An equation relating this ion with constants and other ions
                'value': (sympy.expr | float) Can either be an expression of C at equilibrium of this ion or C itself
                'const_symbol': (str) The constant related to this complex
                'origin': (str) How this complex is related to other ion. Set to 'complex'

        3) (dict) solid_reference
            [key] structure:
                'x' + 'cat' + {cation_id} + 'ani' + {anion_id} + 'sl'
            [value] structure:
                'formula': (str) The formula of this solid (e.g: 'BaSO3')
                'cation': (str) The formula of the cation of this solid
                'anion': (str) The formula of the anion of this solid
                'cation_symbol': (str) The symbol of the cation of this solid, for easier access
                'anion_symbol': (str) The symbol of the anion of this solid, for easier access
                'cat_id': (int) The id of the cation of this solid
                'ani_id': (int)The id of the anion of this solid
                'cation_count': (int) How many cation is in one solid molecule
                'anion_count': (int) How many anion is in one solid molecule
                'eq': (sympy.Eq) An equation relating the ions that made up this solid
                'value': (sympy.expr | float) Can either be an expression of C at equilibrium of this solid or C itself
                'const_symbol': (str) The constant related to this solid
                'origin': (str) How this solid is related to other ion. Set to 'solid'

        4) (dict) init_anion_concentration / (dict) init_cation_concentration
            [key] structure:
                'C0' + 'ani'(anion) / 'cat'(cation) + {id}
            [value] structure:
                'value': (float) Value of C₀ (measured in molarity)
                'id': (int) The cation / anion id of the ion species related to this init_conc
                'species': (str) The ion most characteristic of this species
                'eq': (sympy.Eq) The equation relating init_conc with other ions

        5) (dict) constants
            [key] structure:
                (for anion): 'pKaani' + {anion_id} + 'st' + {conjugated_state}
                (for cation): 'pBcat' + {cation_id} + 'st' + {conjugated_state}
                (for H⁺ and OH⁻) : 'pKw'
            [value] structure:
                'eq': (sympy.Eq) The equation relating this constant and some other ions
                'symbol': Blank for now
                'value': (float) The value of this constant
                'type': (str) pB for 'cation' / pKa for 'anion'
                'id': (int) The anion id / cation id of the ion species
                'cat_id': (int) (complex only) The cation id
                'ani_id': (int) (complex only) The anion id
                'state': (int) Conjugated state k
                'origin': (str) How this constant is related to other ion ['acid-base', 'complex', 'special']
                    If 'origin' is 'special', there is additional info on this constant's formula after the comma
        """
        # Initialize self.constants, self.init_concentration, self.reference
        self.constants.update({'pKw': {'value': '14', 'eq': None, 'symbol': None, 'type': 'cation', 'id': None, 'origin': 'acid-base'}})
        self.cation_reference.update({'h': {'value': None, 'charge': 1, 'id': None, 'formula': 'H+', 'method': None, 'eq': None,
                                 'type': 'cation', 'multiplier': 1, 'state': 0,
                                            'const_symbol': None, 'origin': 'acid-base'}})

        # Add conjugated base anions to self.anion_reference, update self.init_anion_concentration
        for i, anion in enumerate(self.anions.values()):  # anion: matter.Anion object
            # Skip OH- because it is irrelevant ( we already added it at the start)
            if anion.formula == 'OH-':
                continue
            self.init_anion_concentration.update(
                {f'C0ani{i}': {'species': anion.formula, 'id': i, 'value': anion.init.M, 'eq': None}}
            )
            # Add conjugated base form (state = k) of anion to reference
            self.anion_reference.update({
                    f'xani{i}st{anion.conjugate_count}': {'id': i, 'value': None, 'eq': None,
                                'state': anion.conjugate_count,
                                 'charge': anion.charge,
                                'multiplier': anion.multiplier,
                                 'formula': anion.formula,
                                 'origin': 'acid-base',
                                 'method': None, 'const_symbol': f'pKaani{i}st{anion.conjugate_count}'}
            })
            self.logbook[anion.formula] = f'xani{i}st{anion.conjugate_count}'

            # Add conjugated acid forms (state = 0 -> k-1) of anion to reference
            for k in range(anion.conjugate_count):
                # conjugated_form_anion_k is the conjugated form (state = k) of the extreme state = n
                conjugated_form_anion_stk = matter.Anion(datadict=iloc(by=anion.formula, i=k))
                self.anion_reference.update({
                    f'xani{i}st{k}': {'id': i, 'value': None, 'state': k, 'eq': None,
                                      'origin': 'acid-base',
                                       'charge': conjugated_form_anion_stk.charge,
                                      'multiplier': conjugated_form_anion_stk.multiplier,
                                        'formula': conjugated_form_anion_stk.formula,
                                       'const_symbol': f'pKaani{i}st{k}' if k != 0 else None}
                })
                self.logbook[conjugated_form_anion_stk.formula] = f'xani{i}st{k}'
                # Update pKa value
                self.constants.update({f'pKaani{i}st{k+1}':
                                           {'value': f'{anion["pKa"][k+1]}',
                                            'eq': None,
                                            'symbol': None,
                                            'type': 'anion',
                                            'id': i,
                                            'origin': 'acid-base',
                                            'state': k+1}})

            # Count the number of special form of this series (if any)
            special_count = 0 if anion.specials is None else len(anion.specials)
            # Updates all special form. If special_count = 0 (no special form), skips this step entirely.
            for k, special_form_k in enumerate(range(special_count)):
                special_form_k = k + 1
                # special_form_k = k + 1, since numbering for compound starts at 0
                special_form_formula = anion.specials[k].formula
                # Update reference on this special form
                self.anion_reference.update({
                    f'xani{i}spec{special_form_k}': {
                        'id': i, 'value': None, 'eq': None, 'state': special_form_k,
                        'charge': anion.specials[k].charge, 'multiplier': anion.specials[k].multiplier,
                        'formula': anion.specials[k].formula, 'origin': 'special',
                        'const_symbol': f'pKaani{i}spec{special_form_k}',
                    }
                })
                self.logbook[special_form_formula] = f'xani{i}spec{special_form_k}'

                # Also update constant pKspecial
                self.constants.update({f'pKaani{i}spec{special_form_k}':{
                    'value': anion.specials[k].pKspec[special_form_k],
                    'eq': None,
                    'symbol': None,
                    'type': 'anion',
                    'origin': f'special, {anion.specials[k].relation}',
                    'id': i,
                    'state': special_form_k}})

        # Add conjugated acid cations to self.cation_reference, update self.init_cation_concentration
        for i, cation in enumerate(self.cations.values()):  # cation: matter.Cation object
            # Skip H+ because it is irrelevant ( we already added it at the start)
            if cation['formula'] == 'H+':
                continue
            self.init_cation_concentration.update(
                {f'C0cat{i}': {'species': cation.formula, 'id': i, 'value': cation.init.M, 'eq': None}}
            )
            self.cation_reference.update(
                {
                    f'xcat{i}st0': {'id': i, 'value': None, 'eq': None, 'state': 0,
                                    'multiplier': cation.multiplier,
                                    'origin': 'acid-base',
                                 'charge': cation.charge,
                                 'formula': cation.formula,
                                 'method': None, 'const_symbol': None}
                }
            )
            self.logbook[cation['formula']] = f'xcat{i}st0'
            # Add other conjugated form of cations to self.cation_reference
            for k in range(cation['conjugate_count']):
                # Since acid start at state = n, but cation start at state = 0, so we skip k = 0 (alreaedy added state 0)
                k = k + 1
                conjugated_form_cation_stk = matter.Cation(datadict=iloc(by=cation.formula, i=k))
                self.cation_reference.update({
                    f'xcat{i}st{k}': {'id': i, 'value': None, 'state': k, 'eq': None,
                                      'origin': 'acid-base',
                                      'multiplier': conjugated_form_cation_stk.multiplier,
                                       'charge': conjugated_form_cation_stk.charge,
                                       'formula': conjugated_form_cation_stk.formula,
                                       'const_symbol': f'pBcat{i}st{k}' if k != 0 else None}
                })
                self.logbook[conjugated_form_cation_stk['formula']] = f'xcat{i}st{k}'
                # Update pB value, except when k = 0 (there is no pB0)
                if k != 0:
                    self.constants.update({f'pBcat{i}st{k}':
                                               {'value': f'{cation.pKa[k]}',  # TODO: Maybe change this back to pB?
                                                'eq': None,
                                                'symbol': None,
                                                'type': 'cation',
                                                'origin': 'acid-base',
                                                'id': i,
                                                'state': k}})

        # Add complexes to self.complex_referemce
        cation_list = [datadict['formula'] for datadict in self.cation_reference.values()]
        anion_list = [datadict['formula'] for datadict in self.anion_reference.values()]
        cation_anion_pair_list = list(product(cation_list, anion_list))
        for cation_anion_pair in cation_anion_pair_list:
            curr_catani_complex_dict = loc(by=list(cation_anion_pair), chemtype='complex')
            if len(curr_catani_complex_dict) > 0:
                for comp_formula, comp_data in curr_catani_complex_dict.items():
                    complex_cat_id = self.cation_reference[self.symbol_of(comp_data['cation'])]['id']
                    complex_ani_id = self.anion_reference[self.symbol_of(comp_data['ligand'])]['id']
                    comp_symbol = f'xcat{complex_cat_id}ani{complex_ani_id}cp{comp_data["ligand_count"]}'
                    self.complex_reference[comp_symbol] = {
                        'formula': comp_formula,
                        'cation': comp_data['cation'],
                        'ligand': comp_data['ligand'],
                        'cation_symbol': self.symbol_of(comp_data['cation']),
                        'ligand_symbol': self.symbol_of(comp_data['ligand']),
                        'cat_id': complex_cat_id,
                        'ani_id': complex_ani_id,
                        'charge': comp_data['charge'],
                        'cation_count': comp_data['cation_count'],
                        'ligand_count': comp_data['ligand_count'],
                        'state': comp_data['ligand_count'],
                        'eq': None,
                        'value': None,
                        'const_symbol': f'logBcat{complex_cat_id}ani{complex_ani_id}cp{comp_data["ligand_count"]}',
                        'origin': 'complex'
                    }
                    self.logbook[comp_formula] = comp_symbol
                    self.constants[f'logBcat{complex_cat_id}ani{complex_ani_id}cp{comp_data["ligand_count"]}'] = {
                        'value': comp_data['logB'][comp_data['ligand_count']],  # ligand_count is also state i
                        'eq': None,
                        'symbol': None,
                        'type': 'complex',
                        'origin': 'complex',
                        'cat_id': complex_cat_id,
                        'ani_id': complex_ani_id,
                        'state': comp_data['ligand_count']}

        # Add solids to self.solid_reference
        for cation_anion_pair in cation_anion_pair_list:
            curr_solid_dict = loc(by=[cation_anion_pair[0], cation_anion_pair[1]], chemtype='solid')
            if len(curr_solid_dict) > 0:
                solid_cat_id = self.cation_reference[self.symbol_of(curr_solid_dict['cation'])]['id']
                solid_ani_id = self.anion_reference[self.symbol_of(curr_solid_dict['anion'])]['id']
                solid_symbol = f'xcat{solid_cat_id}ani{solid_ani_id}sl'
                self.solid_reference[solid_symbol] = {
                    'formula': curr_solid_dict['formula'],
                    'cation': curr_solid_dict['cation'],
                    'anion': curr_solid_dict['anion'],
                    'cation_symbol': self.symbol_of(curr_solid_dict['cation']),
                    'anion_symbol': self.symbol_of(curr_solid_dict['anion']),
                    'cat_id': solid_cat_id,
                    'ani_id': solid_ani_id,
                    'cation_count': curr_solid_dict['cation_count'],
                    'anion_count': curr_solid_dict['anion_count'],
                    'eq': None,
                    'value': None,
                    'const_symbol': f'pKscat{solid_cat_id}ani{solid_ani_id}sl',
                    'origin': 'solid'
                }
                self.logbook[solid_symbol] = curr_solid_dict['formula']
                self.constants[f'pKscat{solid_cat_id}ani{solid_ani_id}sl'] = {
                    'value': curr_solid_dict['pKs'][1], # TODO: I put a 1 here for now which is prob correct for all solids anyways
                    'eq': None,
                    'symbol': None,
                    'type': 'solid',
                    'origin': 'solid',
                    'cat_id': solid_cat_id,
                    'ani_id': solid_ani_id,
                    'state': 0}

    def _construct_const_eq(self) :
        for const, const_info in self.constants.items():
            # Skip pKw
            if const == 'pKw':
                continue
            # Cation and anion has the same structure
            # Construct the symbol for HA and A-, or A+ and AOH
            if const_info['origin'] == 'acid-base':
                if const_info['type'] == 'anion':
                    conj_base = self.construct_symbol(const_info['id'], const_info['state'], const_info['type'])
                    conj_acid = self.construct_symbol(const_info['id'], const_info['state'] - 1, const_info['type'])
                    # Create this type of equation: Ka = [H+][A-]/[HA]
                    eq_r = f'(h * {conj_base}) / {conj_acid}'
                    eq_l = f'10 ** -{const}'
                if const_info['type'] == 'cation':
                    conj_base = self.construct_symbol(const_info['id'], const_info['state'], const_info['type'])
                    conj_acid = self.construct_symbol(const_info['id'], 0, const_info['type'])
                    # Create this type of equation: 10^-pB2 = [H+]^2 [A(OH)2]/[A2+]
                    eq_r = f'(h ** {const_info["state"]} * {conj_base}) / {conj_acid}'
                    eq_l = f'10 ** -{const}'

            # For ions formed in a special way
            if 'special' in const_info['origin']:
                # eq_r is something like '[h] * [spec] / ([st0]**2)'. We need to convert this into an equation
                eq_l = f'10 ** -{const}'
                eq_r = const_info['origin'].split(',')[-1].strip()
                # This replaces [] with usable symbols
                eq_r = eq_r.replace('[h]', 'h').replace('[spec]', f'xani{const_info["id"]}spec{const_info["state"]}')\
                    .replace('[st0]', f'xani{const_info["id"]}st0').replace('[st1]', f'xani{const_info["id"]}st1')\
                    .replace('[st2]', f'xani{const_info["id"]}st2').replace('[st3]', f'xani{const_info["id"]}st3')

            # For complexes
            if const_info['origin'] == 'complex':
                comp_symbol = f'xcat{const_info["cat_id"]}ani{const_info["ani_id"]}cp{const_info["state"]}'
                comp_cation = self.complex_reference[comp_symbol]['cation_symbol']
                comp_ligand = self.complex_reference[comp_symbol]['ligand_symbol']
                ligand_count = self.complex_reference[comp_symbol]['ligand_count']

                # Create this type of equation: 10^logB2 = [A+][B-]^2/[AB2-]
                eq_r = f'{comp_symbol} / ({comp_cation} * {comp_ligand} ** {ligand_count})'
                eq_l = f'10 ** {const}'

            # For solids
            if const_info['origin'] == 'solid':
                solid_symbol = f'xcat{const_info["cat_id"]}ani{const_info["ani_id"]}sl'
                solid_cation = self.solid_reference[solid_symbol]['cation_symbol']
                solid_anion = self.solid_reference[solid_symbol]['anion_symbol']
                cation_count = self.solid_reference[solid_symbol]['cation_count']
                anion_count = self.solid_reference[solid_symbol]['anion_count']

                # Create this type of equation: 10^-pKs = [A+]^a * [B-]^b
                eq_r = f'({solid_cation} ** {cation_count}) * ({solid_anion} ** {anion_count})'
                eq_l = f'10 ** -{const}'

            # FINALLY, Change string to sympy expression
            eq_r = sympy.parse_expr(eq_r)
            eq_l = sympy.parse_expr(eq_l)
            eq = sympy.Eq(eq_l, eq_r)
            # Update value in self.constants
            self.constants[const]['eq'] = eq

    def _construct_init_M_eq(self):
        # For anion (C0ani)
        for init_M_symbol, init_M_info in self.init_anion_concentration.items():
            # Get a dictionary of all ions related to this C0 along with its multiplier (count of elem)
            dict_kspecies_valmultiplier = self.fetch_related(ani_id=init_M_info['id'],
                                                             dict_return={'anion': 'multiplier', 'complex': 'ligand_count',
                                                                          'solid': 'anion_count'})
            eq_l = f'{init_M_symbol}'
            eq_r = f''
            for ion, multiplier in dict_kspecies_valmultiplier.items():
                # estimate
                # Sample all anions with state = 0 and check pKa1. If it is < -3, skip this anion (estimation)
                # if self.anion_reference.get(ion, {'state': -1})['state'] == 0:
                #     if float(self.constants[f'pKaani{self.anion_reference[ion]["id"]}st1']['value']) < -3:
                #         continue
                # Reminder: An example of eq_r: 1 * xani0st0 + 1 * xani0st1
                # If multiplier = 1, we don't need to include it because 1 * xani0st0 = xani0st0
                if multiplier == 1:
                    eq_r = join(eq=eq_r, new=ion, operator='+')
                elif multiplier == 0:
                    continue
                else:
                    eq_r = join(eq=eq_r, new=multiplier, operator='+')
                    eq_r = join(eq=eq_r, new=ion, operator='*')
            # Testing out new functionality: replace eq_l with the value because why not
            eq_l = f'{init_M_info["value"]}'
            eq_r = sympy.parse_expr(eq_r)
            eq_l = sympy.parse_expr(eq_l)
            eq = sympy.Eq(eq_l, eq_r)
            self.init_anion_concentration[init_M_symbol].update({'eq': eq})

        # For cation (C0cat) (should be the same)
        for init_M_symbol, init_M_info in self.init_cation_concentration.items():
            # Get a dictionary of all ions related to this C0 along with its multiplier (count of elem)
            dict_kspecies_valmultiplier = self.fetch_related(cat_id=init_M_info['id'],
                                                             dict_return={'cation': 'multiplier',
                                                                          'complex': 'cation_count',
                                                                          'solid': 'cation_count'})
            eq_l = f'{init_M_symbol}'
            eq_r = f''
            for ion, multiplier in dict_kspecies_valmultiplier.items():
                # Reminder: An example of eq_r: 1 * xani0st0 + 1 * xani0st1
                # If multiplier = 1, we don't need to include it because 1 * xani0st0 = xani0st0
                if multiplier == 1:
                    eq_r = join(eq=eq_r, new=ion, operator='+')
                elif multiplier == 0:
                    continue
                else:
                    eq_r = join(eq=eq_r, new=multiplier, operator='+')
                    eq_r = join(eq=eq_r, new=ion, operator='*')
            # Testing out new functionality: replace eq_l with the value because why not
            eq_l = f'{init_M_info["value"]}'
            eq_r = sympy.parse_expr(eq_r)
            eq_l = sympy.parse_expr(eq_l)
            eq = sympy.Eq(eq_l, eq_r)
            self.init_cation_concentration[init_M_symbol].update({'eq': eq})

    def _construct_ion_value_eq(self):
        # For anion (state > 0)
        for ani_sym, ani_info in self.anion_reference.items():
            # Equation for state > 0 (related through pK / pB value), special is uses the same mechanism
            if ani_info['state'] > 0 or ani_info['origin'] == 'special':
                # For state > 0, we can get its equation by solving for its value using the constant pKa equation
                const_eq = self.constants[ani_info['const_symbol']]['eq']
                # Sometimes evaluating solveset for special values lead to a complement instead of a finite set
                # So we test if it is a complement because next(iter()) can't be used on a complement for whatever reason
                ani_sym_val = sympy.solveset(const_eq, sympy.Symbol(f'{ani_sym}'))
                if isinstance(ani_sym_val, sympy.Complement):
                    ani_sym_val = ani_sym_val.as_relational(sympy.Symbol(ani_sym)) # Turn into And object
                    ani_sym_val = ani_sym_val.args # Turn And object into tuple for extraction
                    # Iterating over tuple and take only the equality
                    for val in ani_sym_val:
                        if isinstance(val, Equality):
                            ani_sym_val = val
                            break
                    ani_sym_val = ani_sym_val.args # Turn Equality object into tuple
                    ani_sym_val = ani_sym_val[1] # Item on the right is the expr we are looking for, the left is symbol
                    self.anion_reference[ani_sym]['value'] = ani_sym_val
                else:
                    self.anion_reference[ani_sym]['value'] = next(iter(ani_sym_val))
                eq_l = sympy.parse_expr(ani_sym)
                eq_r = self.anion_reference[ani_sym]['value']
                self.anion_reference[ani_sym]['eq'] = sympy.Eq(eq_l, eq_r)

        # For cation (state > 0)
        for cat_sym, cat_info in self.cation_reference.items():
            # Equation for state > 0 (related through pK / pB value)
            # Skip H+
            if cat_sym == 'h':
                continue
            if cat_info['state'] > 0:
                # For state > 0, we can get its equation by solving for its value using the constant pB equation
                const_eq = self.constants[cat_info['const_symbol']]['eq']
                self.cation_reference[cat_sym]['value'] = next(iter(sympy.solveset(const_eq,
                                                                                  sympy.Symbol(f'{cat_sym}'))))
                eq_l = sympy.parse_expr(cat_sym)
                eq_r = self.cation_reference[cat_sym]['value']
                self.cation_reference[cat_sym]['eq'] = sympy.Eq(eq_l, eq_r)

        # For complex
        for cp_sym, cp_info in self.complex_reference.items():
            subs_dict = self.fetch(item='reference', chemtype=['anion', 'cation'], dict_return='value')
            related_const_eq = self.constants[cp_info['const_symbol']]['eq']
            done = False
            # Exhaustive substitution until there are no state k in equation
            while not done:
                t = related_const_eq.subs(subs_dict)
                if related_const_eq == t:
                    done = True
                else:
                    related_const_eq = t
            self.complex_reference[cp_sym]['value'] = next(iter(sympy.solveset(related_const_eq,
                                                                               sympy.Symbol(f'{cp_sym}'))))
            eq_l = sympy.parse_expr(cp_sym)
            eq_r = self.complex_reference[cp_sym]['value']
            cp_eq = sympy.Eq(eq_l, eq_r)
            self.complex_reference[cp_sym]['eq'] = cp_eq

        # For anion (state = 0)
        # We have to loop twice because the equation for state 0 requires other states value calculated beforehand
        for ani_sym, ani_info in self.anion_reference.items():
            # Equation for state 0, we can get by solving for its value using C0 equation, then exhaustively substitute
            if ani_info['state'] == 0 and ani_info['origin'] == 'acid-base':
                # Get substitute dict by setting sol_id = -1 for fetch
                subs_dict = self.fetch(item='reference', chemtype=['anion', 'complex'], dict_return='value')
                init_M_eq = self.init_anion_concentration[f'C0ani{ani_info["id"]}']['eq']
                done = False
                # Exhaustive substitution until there are no state k in equation
                while not done:
                    t = init_M_eq.subs(subs_dict)
                    if init_M_eq == t:
                        done = True
                    else:
                        init_M_eq = t
                # Update value
                self.equations.update({ani_sym: init_M_eq})
                # Solve for state 0 and update reference. However, there may be more than 1 possible solution, one
                # larger than 0 and one smaller than 0. Obviously we are looking for the solution where state0 > 0.
                # So we call extension.SympyExtension.ensure_pos() to find the positive solution
                sol = list()
                try:
                    for possible_sol in list(sympy.solveset(init_M_eq, sympy.Symbol(ani_sym))):
                        if extension.SympyExtension.ensure_pos(possible_sol):
                            sol.append(possible_sol)
                    assert len(sol) == 1
                    self.anion_reference[ani_sym]['value'] = sol[0]
                    self.anion_reference[ani_sym]['eq'] = sympy.Eq(sympy.Symbol(ani_sym),
                                                                  self.anion_reference[ani_sym]['value'])
                except Exception as e:
                    warnings.warn('Exhaustive substitution not possible')

        # For cation (state = 0)
        for cat_sym, cat_info in self.cation_reference.items():
            # Equation for state 0, we can get by solving for its value using C0 equation, then exhaustively substitute
            # Skip H+
            if cat_sym == 'h':
                continue
            if cat_info['state'] == 0:
                # Get substitute dict by setting sol_id = -1 for list_all_from_id
                subs_dict = self.fetch(item='reference', chemtype=['cation', 'complex'], dict_return='value')
                init_M_eq = self.init_cation_concentration[f'C0cat{cat_info["id"]}']['eq']
                done = False
                # Exhaustive substitution until there are no state k in equation
                while not done:
                    t = init_M_eq.subs(subs_dict)
                    if init_M_eq == t:
                        done = True
                    else:
                        init_M_eq = t
                # Update value
                self.equations.update({cat_sym: init_M_eq})
                try:
                    self.cation_reference[cat_sym]['value'] = next(iter(sympy.solveset(init_M_eq,
                                                                                      sympy.Symbol(cat_sym))))
                    self.cation_reference[cat_sym]['eq'] = sympy.Eq(sympy.Symbol(cat_sym),
                                                                   self.cation_reference[cat_sym]['value'])
                except Exception as e:
                    warnings.warn('Exhaustive substitution not possible')

    def _construct_solid_eq(self):
        subs_dict_1 = {key: val['value'] for key, val in self.anion_reference.items() if not (val['state'] == 0 and
                                                                                              val[
                                                                                                  'origin'] == 'acid-base')}
        subs_dict_2 = {key: val['value'] for key, val in self.cation_reference.items() if not (val['state'] == 0 and
                                                                                               val[
                                                                                                   'origin'] == 'acid-base')}
        subs_dict = {**subs_dict_1, **subs_dict_2}
        for solid_symbol, solid_info in self.solid_reference.items():
            # FULL SUBSTITUION, STATE 0 EXCLUDED
            done = False
            eq = self.constants[solid_info['const_symbol']]['eq']
            while not done:
                t = eq.subs(subs_dict)
                if eq == t:
                    done = True
                else:
                    eq = t
            self.equations[solid_symbol] = eq

    def _construct_ion_conservation_eq(self, substitute = True):
        """
        Will substitute values by default
        :param reference: reference
        :return:
        """
        # eq_l is for positive ion, left empty at first
        eq_l = ''
        # eq_r is for negative ion, start with OH-
        eq_r = '(1e-14)/h'
        # We iterate all ions and build up the ion conservation equation
        dict_kions_valcharges = self.fetch(item='reference', chemtype=['cation', 'anion', 'complex'], dict_return='charge')
        for ion_sym, ion_charge in dict_kions_valcharges.items():
            if ion_charge == 0:
                continue
            if ion_charge > 0:
                eq_l = join(eq=eq_l, new=ion_charge, operator='+')
                eq_l = join(eq=eq_l, new=ion_sym, operator='*')
            if ion_charge < 0:
                eq_r = join(eq=eq_r, new=abs(ion_charge), operator='+')
                eq_r = join(eq=eq_r, new=ion_sym, operator='*')
        eq_l = sympy.parse_expr(eq_l)
        eq_r = sympy.parse_expr(eq_r)
        eq = sympy.Eq(eq_l, eq_r)
        self.ion_eq_tier_1 = eq

        # FULL SUBSTITUION, STATE 0 EXCLUDED
        subs_dict_1 = {key: val['value'] for key, val in self.anion_reference.items() if not (val['state'] == 0 and
                                                                                        val['origin'] == 'acid-base')}
        subs_dict_2 = {key: val['value'] for key, val in self.cation_reference.items() if not (val['state'] == 0 and
                                                                                        val['origin'] == 'acid-base')}
        subs_dict_3 = {key: val['value'] for key, val in self.complex_reference.items()}
        subs_dict = {**subs_dict_1, **subs_dict_2, **subs_dict_3}
        done = False
        eq = self.ion_eq_tier_1
        while not done:
            t = eq.subs(subs_dict)
            if eq == t:
                done = True
            else:
                eq = t
        self.equations.update({'h': eq})

    def _estimate_state_0(self, reference) -> dict:
        for state_symbol, state_info in reference['x'].items():
            # LEGACY CODE BUT I HAVENT FIXED IT YET
            if state_info.get('state', -1) == 0:
                reference['x'][state_symbol]['estimation'] = \
                    reference['init_M'][f'ΣCsol{state_info["sol_id"]}']['value']
            # END
            # We don't estimate molarity of H+
            if state_symbol == 'h':
                continue
            estimation = state_info['value']
            # Get estimation of molarity when pH = 1 (C = 0.1M)
            _sub_dict = dict(self._val_dict)
            _sub_dict.update({'h': 0.1})
            estimate_pH_1 = sympy.N(estimation.subs(_sub_dict), 1)
            # Get estimation of molarity when pH = 14 (C = 1e-14M)
            _sub_dict = dict(self._val_dict)
            _sub_dict.update({'h': 1e-14})
            estimate_pH_14 = sympy.N(estimation.subs(_sub_dict), 1)
            # If pH value at two extremities are near each other, set its value and eq equal to its estimation at pH=7
            if len(estimate_pH_14.free_symbols) == 0:
                if estimate_pH_1 == estimate_pH_14:
                    _sub_dict = dict(self._val_dict)
                    _sub_dict.update({'h': 1e-7})
                    estimate_pH_7 = sympy.N(estimation.subs(_sub_dict), 5)
                    reference['x'][state_symbol]['value'] = estimate_pH_7
                    reference['x'][state_symbol]['eq'] = sympy.Eq(sympy.parse_expr(state_symbol), estimate_pH_7)
        return reference

    def pH_solve_wolfram(self):
        # Initialize Wolfram solver
        wolfram_solver = WolframSolver()

        # Declare all substituable symbols (constants)
        dict_kconstants_valvalue = self.fetch(item='constant', sol_id=-1, dict_return='value')
        for sym, value in dict_kconstants_valvalue.items():
            wolfram_solver.write_expression(sym, value)

        # Declare all equations, different depending on which solve type it is
        for sym, eq in self.equations.items():
            wolfram_solver.write_equation(eq)
        # Tell matlab solver which symbol to solve for
        syms_to_solve = list(self.equations.keys())

        index_of_h = syms_to_solve.index('h')
        syms_to_solve[0], syms_to_solve[index_of_h] = syms_to_solve[index_of_h], syms_to_solve[0]
        wolfram_solver.syms_to_solve = syms_to_solve

        # Solve!
        pH = wolfram_solver.compute_eq()
        return pH

    def pH_solve_gekko(self):
        # Initialize Wolfram solver
        gekko_solver = GEKKOSolver(strict_check=self.options.STRICT_CHECK)

        # Get all constants
        dict_constant_symbols = self.fetch(item='constant', sol_id=-1, dict_return='value')

        # Get all variables
        list_variable_symbols = list(self.equations.keys())

        # Get all equations in list form
        list_eqns = [eqn for eqn in self.equations.values()]

        # Solve!
        pH = gekko_solver.compute_eq(eqns=list_eqns, variables=list_variable_symbols, constants=dict_constant_symbols)
        return pH

    def symbol_of(self, formula: str):
        """
        Returns the symbol (str) given an ion formula
        E.g: symbol_of('Na+') -> 'xcat0st0'    symbol_of('PO43-') -> 'xani0st3'
        :param (str) formula:
        :return:
        """
        return self.logbook[formula]  # I don't use get to catch exceptions.

    def construct_symbol(self, sol_id: int, state: int, chemtype: str) -> str:
        """
        Quickly construct a symbol from given parameters.
        E.g: construct_symbol(0, 1, 'anion') -> 'xani0st1'
        :param sol_id:
        :param state:
        :param type:
        :return:
        """
        if chemtype == 'anion':
            return f'xani{sol_id}st{state}'
        if chemtype == 'cation':
            return f'xcat{sol_id}st{state}'

    def fetch(self, item: str | list = ['reference', 'init_concentration', 'constant'],
              sol_id: int = -1, origin: str = '',
              chemtype: str| list = ['cation', 'anion'], dict_return: str | dict = None,
              exclude_sol_id: int | list = [], exclude_state: int | list = [],
              skip_empty: bool = True) -> list | dict:
        """
        If dict_return is a list of Nones, return a list of all items requested given its sol_id and chemtype
        E.g: fetch(sol_id=0, chemtype=['anion']) -> ['xani0st0', 'xani0st1', 'C0ani0', 'pKaani0st1']
        E.g: fetch(item=['reference'], sol_id=0, chemtype=['anion']) -> ['xani0st0', 'xani0st1']

        If dict_return is a valid key,return a dictionary of all items requested given its sol_id and chemtype,
        its value will be set to what dict_return respectively
        E.g: fetch(sol_id=0, chemtype='anion', dict_return=['reference': 'multiplier', 'init_concentration': 'multiplier',
        'constant': 'multiplier']) -> {'xani0st0': 1, 'xani0st1': 1}

        If sol_id is equal to -1, this will select all ion regardless of id
        If chemtype is equal to '', this will select all chemtype

        :param (list) item: Which item will be returned. Must be in ['reference', 'init_concentration', 'constant'].
            Select all by default
            Alternatively, if there is only 1 item, you can also pass only 1 string. E.g: 'reference'
        :param (int) sol_id: Id of the species. Set equal to -1 to select all. Select all by default
        :param (str) chemtype: A list of criteria for chemtype. Select all by default
            'anion' for anion and 'cation' for cation. Set to ['cation', 'anion'] to select all ion
        :param (str) dict_return: The value to set for each item (key)
             E.g: {'reference': 'value', 'init_concentration': 'value', 'constant': 'value'}
             Alternatively, if there is only 1 value, you can pass in a string. E.g: 'value'
             By default, set to None (returns list)
             All the keys of self.reference are valid dict_return arguments (see below)
             [formula, id, charge, state, eq, value, multiplier, const_symbol, method]
        :param (list) exclude_sol_id: A list of sol_id to ignore while fetching. Empty by default.
        :param (bool) skip_empty: For dictionary fetching only. If a respected item's value is empty, skip it entirely.
            True by default
        :return: list | dict
        """
        # Input checking, handling
        # If item is a string, convert it to a list of 1 string
        if isinstance(item, str):
            item = [item]
        # If exclude_sol_id is an integer, convert it to a list of 1 integer
        if isinstance(exclude_sol_id, int):
            exclude_sol_id = [exclude_sol_id]
        # If exclude_state is an integer, convert it to a list of 1 integer
        if isinstance(exclude_state, int):
            exclude_state = [exclude_state]
        # If dict_return is a string, convert it to a dictionary of keys = items and value = that string
        if isinstance(dict_return, str):
            dict_return = {key: dict_return for key in item}
        # If chemtype is a string, convert it to a list of 1 string
        if isinstance(chemtype, str):
            chemtype = [chemtype]

        """List return mode"""
        if dict_return is None:
            return_list = []
            # Returns reference?
            if 'reference' in item:
                # For anion
                if 'anion' in chemtype:
                    for ani_sym, ani_info in self.anion_reference.items():
                        if (sol_id == ani_info['id'] or sol_id == -1) and (ani_info['id'] not in exclude_sol_id) and \
                                (ani_info['state'] not in exclude_state):
                            if origin == ani_info['origin'] or origin == '':
                                return_list.append(ani_sym)
                # For cation
                if 'cation' in chemtype:
                    for cat_sym, cat_info in self.cation_reference.items():
                        if (sol_id == cat_info['id'] or sol_id == -1) and (cat_info['id'] not in exclude_sol_id) and \
                                (cat_info['state'] not in exclude_state):
                            if origin == cat_info['origin'] or origin == '':
                                return_list.append(cat_sym)
                # For complex
                if 'complex' in chemtype:
                    for cp_sym, cp_info in self.complex_reference.items():
                        return_list.append(cp_sym)
            # Returns constants?
            if 'constant' in item:
                for const_sym in self.constants.keys():
                    if self.constants[const_sym]['id'] not in exclude_sol_id:
                        if origin in self.constants[const_sym]['origin'] or origin == '':
                            return_list.append(const_sym)
            # Returns init_concentration
            if 'init_concentration' in item:
                # For anion
                if 'anion' in chemtype:
                    for init_M, init_M_info in self.init_anion_concentration.items():
                        if (init_M_info['id'] not in exclude_sol_id) and (sol_id == init_M_info['id'] or sol_id == -1):
                            return_list.append(init_M)
                # For cation
                if 'cation' in chemtype:
                    for init_M, init_M_info in self.init_cation_concentration.items():
                        if (init_M_info['id'] not in exclude_sol_id) and (sol_id == init_M_info['id'] or sol_id == -1):
                            return_list.append(init_M)
            return return_list

        """Dict return mode"""
        if isinstance(dict_return, dict):
            return_dict = dict()
            # Returns reference?
            if 'reference' in item:
                # For anion
                if 'anion' in chemtype:
                    for ani_sym, ani_info in self.anion_reference.items():
                        if ani_info[dict_return['reference']] is None and skip_empty:
                            continue
                        if (sol_id == ani_info['id'] or sol_id == -1) and (ani_info['id'] not in exclude_sol_id) and \
                                (ani_info['state'] not in exclude_state):
                            if origin == ani_info['origin'] or origin == '':
                                return_dict.update({ani_sym: ani_info[dict_return['reference']]})
                # For cation
                if 'cation' in chemtype:
                    for cat_sym, cat_info in self.cation_reference.items():
                        if cat_info[dict_return['reference']] is None and skip_empty:
                            continue
                        if (sol_id == cat_info['id'] or sol_id == -1) and (cat_info['id'] not in exclude_sol_id) and \
                                (cat_info['state'] not in exclude_state):
                            if origin == cat_info['origin'] or origin == '':
                                return_dict.update({cat_sym: cat_info[dict_return['reference']]})
                # For complex
                if 'complex' in chemtype:
                    for cp_sym, cp_info in self.complex_reference.items():
                        if cp_info.get(dict_return['reference'], None) is None and skip_empty:
                            continue
                        return_dict[cp_sym] = cp_info[dict_return['reference']]
            # Returns constants?
            if 'constant' in item:
                for const_sym, const_info in self.constants.items():
                    if const_info[dict_return['constant']] is None and skip_empty:
                        continue
                    if self.constants[const_sym].get('id', None) not in exclude_sol_id:
                        if origin in const_info['origin'] or origin == '':
                            return_dict.update({const_sym: const_info[dict_return['constant']]})
            # Returns init_concentration
            if 'init_concentration' in item:
                # For anion
                if 'anion' in chemtype:
                    for init_M, init_M_info in self.init_anion_concentration.items():
                        if init_M_info[dict_return['init_concentration']] is None and skip_empty:
                            continue
                        if (init_M_info['id'] not in exclude_sol_id) and (sol_id == init_M_info['id'] or sol_id == -1):
                            return_dict.update({init_M: init_M_info[dict_return['init_concentration']]})
                # For cation
                if 'cation' in chemtype:
                    for init_M, init_M_info in self.init_cation_concentration.items():
                        if init_M_info[dict_return['init_concentration']] is None and skip_empty:
                            continue
                        if (init_M_info['id'] not in exclude_sol_id) and (sol_id == init_M_info['id'] or sol_id == -1):
                            return_dict.update({init_M: init_M_info[dict_return['init_concentration']]})
            return return_dict
        raise AttributeError('Something went wrong while fetching.')

    def fetch_related(self, ani_id=None, cat_id=None, dict_return=None) -> dict | list:
        """
        Find all ions that have a given ani_id or cat_id
        :param ani_id:
        :param cat_id:
        :return:
        """
        if dict_return is None:
            dict_return = {'cation': None, 'anion': None, 'complex': None,
                           'solid': None}
        return_list = list()
        return_dict = dict()
        if cat_id is not None:
            for cat_sym, cat_info in self.cation_reference.items():
                if cat_info['id'] == cat_id:
                    if dict_return['cation'] is None:
                        return_list.append(cat_sym)
                    else:
                        return_dict[cat_sym] = cat_info[dict_return['cation']]
        if ani_id is not None:
            for ani_sym, ani_info in self.anion_reference.items():
                if ani_info['id'] == ani_id:
                    if dict_return['anion'] is None:
                        return_list.append(ani_sym)
                    else:
                        return_dict[ani_sym] = ani_info[dict_return['anion']]
        for cp_sym, cp_info in self.complex_reference.items():
            if cp_info['cat_id'] == cat_id or cp_info['ani_id'] == ani_id:
                if dict_return.get('complex', None) is None:
                    return_list.append(cp_sym)
                else:
                    return_dict[cp_sym] = cp_info[dict_return['complex']]
        for sl_sym, sl_info in self.solid_reference.items():
            if sl_info['cat_id'] == cat_id or sl_info['ani_id'] == ani_id:
                if dict_return.get('solid', None) is None:
                    return_list.append(sl_sym)
                else:
                    return_dict[sl_sym] = sl_info[dict_return['solid']]
        if len(return_dict) > 0:
            return return_dict
        return return_list

    def _update_reference_at_eq(self):
        for symbol, value in self.eq_molarity_dict.items():
            if self.anion_reference.get(symbol, None) is not None:
                self.anion_reference[symbol]['value'] = value
            if self.cation_reference.get(symbol, None) is not None:
                self.cation_reference[symbol]['value'] = value
            if self.complex_reference.get(symbol, None) is not None:
                self.complex_reference[symbol]['value'] = value
        subs_dict = {**self.eq_molarity_dict, **self.fetch(['init_concentration', 'constant'], dict_return='value')}
        for ani_sym, infodict in self.anion_reference.items():
            try:
                self.anion_reference[ani_sym]['value'] = sympy.N(self.anion_reference[ani_sym]['value'].subs(subs_dict))
            except AttributeError as e:
                continue
        for cat_sym, infodict in self.cation_reference.items():
            try:
                self.cation_reference[cat_sym]['value'] = sympy.N(self.cation_reference[cat_sym]['value'].subs(subs_dict))
            except AttributeError as e:
                continue
        for cp_sym, infodict in self.complex_reference.items():
            try:
                self.complex_reference[cp_sym]['value'] = sympy.N(self.complex_reference[cp_sym]['value'].subs(subs_dict))
            except AttributeError as e:
                continue

def join(eq: str, new: str, operator: str = '+') -> str:
    """
    Joins an equation (string) to a new section with an operator, adding spaces where needed for visualization.
    E.g: join('h * x1', 'x2', '+') -> 'h * x1 + x2'
         join('', 'h') -> 'h'
    :param (str) eq: The string to add to
    :param (str) new: The new part to be added to the equation
    :param (str) operator: Mathematical operator ('+', '-', '*', '**', '/')
    :return: (str)
    """
    eq = eq.strip()
    if eq == '' or eq is None:
        return str(new)
    eq = eq + f' {operator} ' + str(new).strip()
    return eq

if __name__ == '__main__':
    """
    hcl = matter.Compound(cation='X+', anion='Cl-', cation_count=1, anion_count=1)
    BaOH2 = matter.Compound(cation='Ba2+', anion='OH-', cation_count=1, anion_count=2)
    CH3COOH = matter.Compound(formula='CH3COOH', cation_count=1, anion_count=1)
    NaOH = matter.Compound(cation='Na+', anion='OH-', cation_count=1, anion_count=1)
    H2SO4 = matter.Compound(formula='H2SO4', cation_count=2, anion_count=1)
    HClO4 = matter.Compound(formula='HClO4', cation_count=1, anion_count=1)
    NaCN = matter.Compound(cation='Na+', anion='CN-', cation_count=1, anion_count=1)
    NaHS = matter.Compound(cation='Na+', anion='HSO3-', cation_count=1, anion_count=1)
    NH3 = matter.Compound(formula='NH3', cation_count=0, anion_count=1)"""
    # BaCl2 = matter.Compound(cation='Ba2+', anion='Cl-', cation_count=1, anion_count=2)
    # Na2CO3 = matter.Compound(cation='Na+', anion='CO32-', cation_count=2, anion_count=1)
    # solA = matter.Solution(volume=1)
    # solA.add_compound(BaCl2, init_M=0.2)
    # solA.add_compound(Na2CO3, init_M=0.2)

    # Cr2SO43 = matter.Compound(cation='Cr3+', anion='SO42-', cation_count=2, anion_count=3)
    # HCl = matter.Compound(cation='H+', anion='Cl-', cation_count= 1, anion_count=1)
    # solA = matter.Solution(volume=1)
    # solA.add_compound(HCl, init_M=1)
    # ana = Analyser(solA)
    # ana.analyze_eq(solver='wolfram')
    # pH = ana.get_pH()