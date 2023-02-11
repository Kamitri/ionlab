from enum import StrEnum
import os
import pathlib

class Solver(StrEnum):
    GEKKO = 'gekko'
    # More? There is more, but I don't want to add them in yet.


class Paths:
    APP_NAME = 'Equion'
    WOLFRAM_KERNEL = None
    PATH_DB_XLSX = f'{pathlib.Path(os.path.dirname(os.path.abspath(__file__))).joinpath("IonicEquilibrium").joinpath("db")}'
    PATH_DB_SQL = f'{pathlib.Path(os.path.dirname(os.path.abspath(__file__))).joinpath("IonicEquilibrium").joinpath("db")}'
    FILENAME_DB_XLSX = 'chemicals.xlsx'
    FILENAME_DB_SQL = 'chemicals.db'


class IonicEquilibriumOptions:
    PRECISION = 15
    TABLES = ('acid', 'base', 'cation', 'anion', 'solid', 'complex')
    NUMERIC_DATA = ('pKa', 'pKs', 'logB')
    SPECIES_NAME = 'species'
    SHORTHAND = {'f': 'formula',
                 'an': 'anion', 'ani': 'anion',
                 'cat': 'cat',
                 'an_c': 'anion_count', 'ac': 'anion_count', 'an_count': 'anion_count', 'anion_c': 'anion_count',
                 'cat_c': 'cation_count', 'cc': 'cation_count', 'cat_count': 'cation_count', 'cation_c': 'cation_count'}
    STRICT_CHECK = True