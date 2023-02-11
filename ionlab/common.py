"""
Handy functions and quick access to a few constant to make live ever just so slightly easier
"""
import pandas
import pandas as pd
from fractions import Fraction

PERIODIC_TABLE = {
     1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', '·': '', 12: 'Mg',
     13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr',
     25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br',
     36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd',
     47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La',
     58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er',
     69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au',
     80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac',
     90: 'Th', 91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
     101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds',
     111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'
}

ACCEPTED_CATIONS = ['Ag', 'Al', 'Ba', 'X']


class AnalysisFailedError(Exception):
    def __init__(self, e):
        super().__init__()
        self.e = e

    def __repr__(self):
        return str(self.e)


def is_user_defined(element_symbol) -> bool:
    if element_symbol in PERIODIC_TABLE.values():
        return False
    return True


def parse_compound_state(state: str) -> list:
    """
    :param state: State as defined in db.xlsx
    :return: List of dictionary containing mapped information
    """
    info = []
    all_states = state.splitlines()
    for single_state in all_states:
        _info = list(map(lambda x: x.strip(), single_state.split(',')))
        info.append({
            'state': _info[0],
            f'{_info[1].split("=")[0].strip()}': float(_info[1].split("=")[1].strip()),
            'method': _info[2]
        })
    return info


def identify_state(state: str | int, compound) -> int:
    if isinstance(state, int):
        if 0 <= state <= 4:
            return state
        else:
            raise AttributeError('State must be between 0 and 4.')
    if isinstance(state, str):
        if state == compound.formula:
            return 0
        for i, single_state in enumerate(compound.states):
            if state == single_state['state']:
                return i + 1
    raise Exception('WOT HAPPENED HERE')


def get_charge(compound, state: str | int) -> int:
    state = identify_state(state, compound)
    return compound.charge - state


def get_complex_charge(cation_charge: int, anion_charge: int, anion_count: int):
    return cation_charge + anion_charge * anion_count


class Data:
    def __init__(self, filename: str):
        self.acid = pandas.read_excel(filename, sheet_name='Acid')
        self.salt = pandas.read_excel(filename, sheet_name='Salt')
        self.complex = pandas.read_excel(filename, sheet_name='Complex')

    def locate_acid(self, search_tag: str) -> pd.Series:
        search_tag = search_tag[0].upper() + search_tag[1:]
        result = self.acid.loc[(self.acid['Name'] == search_tag) | (self.acid['Formula'] == search_tag)]
        result = result.squeeze(axis=0)
        try:
            result.values[0]
        except IndexError:
            raise FileNotFoundError('Compound doesn\'t exist.')
        return result

    def locate_metal_cation(self, search_tag: str) -> pd.Series:
        search_tag = search_tag[0].upper() + search_tag[1:]
        result = self.salt.loc[(self.salt['Formula'] == search_tag)]
        result = result.squeeze(axis=0)
        try:
            result.values[0]
        except IndexError:
            raise FileNotFoundError('Compound doesn\'t exist.')
        return result

    def locate_salt(self, cation: str, anion: str):
        """

        :param search_tag: only supports formula-searching. E.g: Al2(SO4)3
        :return:
        """
        # GET CATION INFO
        _cation = self.salt.loc[self.salt['Formula'] == cation]
        if _cation.empty:
            _cation = self.acid.loc[self.acid['Formula'] == cation]
        if _cation.empty:
            _cation = self.acid.iloc[self.acid['State'].str.contains(cation).tolist().index(True)]
        cation = _cation.squeeze(axis=0)

        # GET ANION INFO
        _anion = self.acid.iloc[self.acid['State'].str.contains(anion).tolist().index(True)]
        anion_state = _anion.State.splitlines()
        anion_state = list(map(lambda x: x.split(',')[0], anion_state)).index(anion) + 1
        if _anion.empty:
            _anion = self.acid.loc[self.acid['Formula'] == anion]
        if _anion.empty:
            _anion = self.salt.loc[self.salt['Formula'] == anion]
        anion = _anion.squeeze(axis=0)

        if anion.empty or cation.empty:
            raise FileNotFoundError('Compound doesn\'t exist.')

        # COUNT CONTRIBUTION OF EACH ION
        cation_charge = cation.Charge
        anion_charge = anion.Charge - anion_state
        frac = Fraction(abs(cation_charge), abs(anion_charge))
        cation_count = frac.denominator
        anion_count = frac.numerator

        return {'cation': {'count': cation_count, 'info': cation, 'charge': cation_charge},
                'anion': {'count': anion_count, 'info': anion, 'charge': anion_charge}}

    def locate_complex(self, pair: list):
        result = self.complex.loc[(self.complex['Metal cation'] == pair[0]) & (self.complex['Ligand'] == pair[1])]
        result = result.squeeze(axis=0)
        try:
            result.values[0]
        except IndexError:
            raise FileNotFoundError('Compound doesn\'t exist.')
        result['State'] = parse_compound_state(result['State'])
        result['Charge'] = list(map(int, result['Charge'].split(',')))
        return result

    def get_complex_list(self) -> list:
        _complex = pd.DataFrame(self.complex.drop('State', axis=1))
        _list = []
        for pair in _complex.iterrows():
            _list.append([pair[1]['Metal cation'], pair[1]['Ligand']])
        return _list


DB_PATH = 'db.xlsx'
db = Data(DB_PATH)

if __name__ == '__main__':
    db.get_complex_list()