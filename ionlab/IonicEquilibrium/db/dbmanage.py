import pandas as pd
import sqlite3
import os
import sys
import traceback
import argparse
import warnings
from pprint import pprint
from typing import Union, List
import ionlab.conf as conf


def make_db(path_xlsx: str = conf.Paths.PATH_DB_XLSX, filename_xlsx: str = conf.Paths.FILENAME_DB_XLSX,
            path_sql: str = conf.Paths.PATH_DB_SQL, filename_sql: str = conf.Paths.FILENAME_DB_SQL) -> None:
    """
    Create a SQL database of all chemicals' information, converted from a .xlsx file.
    :param (str) path_xlsx: The path to the .xlsx file used to create the SQL database.
    :param (str) filename_xlsx: The .xlsx file's name.
    :param (str) path_sql: The path to the output .xlsx file.
    :param (str) filename_sql: The SQL file's name.
    :return: None
    """
    full_path_xlsx = os.path.join(path_xlsx, filename_xlsx)
    full_path_sql = os.path.join(path_sql, filename_sql)
    if not os.path.exists(full_path_sql):
        with open(full_path_sql, 'w') as file:
            pass
    conn = sqlite3.connect(full_path_sql)
    # Create the six tables in order: acid, base, cation, anion, solid, complex.
    conn.execute(
        '''CREATE TABLE IF NOT EXISTS acid
        (formula           TEXT,
        name               TEXT,
        pKa                TEXT,
        species            TEXT,
        conjugate_count    INTEGER,
        charge             INTEGER,
        multiplier         INTEGER
        );''')
    conn.execute(
        '''CREATE TABLE IF NOT EXISTS base
        (formula           TEXT,
        name               TEXT,
        pKa                TEXT,
        species            TEXT,
        conjugate_count    INTEGER,
        charge             INTEGER,
        multiplier         INTEGER
        );''')
    conn.execute(
        '''CREATE TABLE IF NOT EXISTS cation
        (formula           TEXT,
        name               TEXT,
        pKa                TEXT,
        species            TEXT,
        conjugate_count    INTEGER,
        charge             INTEGER,
        multiplier         INTEGER
        );''')
    conn.execute(
        '''CREATE TABLE IF NOT EXISTS anion
        (formula           TEXT,
        name               TEXT,
        pKa                TEXT,
        species            TEXT,
        conjugate_count    INTEGER,
        charge             INTEGER,
        proton_count       INTEGER,
        multiplier         INTEGER,
        special            TEXT
        );''')
    conn.execute(
        '''CREATE TABLE IF NOT EXISTS solid
        (formula           TEXT,
        cation             TEXT,
        anion              TEXT,
        pKs                TEXT,
        cation_count       INTEGER,
        cation_charge      INTEGER,
        anion_count        INTEGER,
        anion_charge       INTEGER
        );''')
    conn.execute(
        '''CREATE TABLE IF NOT EXISTS complex
        (formula           TEXT,
        cation             TEXT,
        ligand             TEXT,
        logB               TEXT,
        species            TEXT,
        cation_count       INTEGER,
        ligand_count       INTEGER,
        species_count      INTEGER,
        charge             INT,
        tier               INT
        );''')

    # Import data from xlsx to sql
    sheets = ['acid', 'base', 'cation', 'anion', 'solid', 'complex']
    for sheet in sheets:
        pd.read_excel(full_path_xlsx, sheet_name=sheet).to_sql(sheet, conn, index=False, if_exists='replace')
    conn.commit()
    conn.close()


def clean_data(datadict: dict) -> dict:
    for data_key, data_val in datadict.items():
        if isinstance(data_val, int) or data_val is None:  # int is fine
            continue
        # Correctly format an array of value.
        if data_key in conf.IonicEquilibriumOptions.NUMERIC_DATA:  # (str) '-3, 1.99' -> {1: -3, 2: 1.99}
            datadict[data_key] = {i+1: x for i, x in enumerate(map(lambda x: float(x.strip()),
                                                                   str(data_val).split(',')))}
        elif isinstance(data_val, float) and data_val.is_integer():  # Converts 1.0 to 1
            datadict[data_key] = int(data_val)
        elif data_key == 'special' and data_val is not None:  # special has different spliting method
            datadict[data_key] = list(map(lambda x: x.strip(), data_val.split(';')))
            for i, special_comp in enumerate(datadict[data_key]):
                special_datalist = list(map(lambda x: x.strip(), datadict[data_key][i].split(',')))
                datadict[data_key][i] = {'formula': special_datalist[0],
                                         'pKspec': {i+1: float(special_datalist[1])},
                                         'relation': special_datalist[2]}
        elif isinstance(data_val, str) and ',' in data_val:  # (str) 'H2SO4, HSO4-, SO42-' -> ['H2SO4', 'HSO4-', 'SO42-']
            datadict[data_key] = list(map(lambda x: x.strip(), data_val.split(',')))
    return datadict


def loc(by: Union[str, list], chemtype: str = None, clean: bool = True) -> dict:
    """
    Locates a chemical within a database given its formula and type (optional) and return its info in a dictionary.
    For chemtype = 'solid' or 'complex', list searching is supported by inputting a list of one cation and one anion.
    If list searching is selected and chemtype wasn't given, the function will look for solid automatically.
    :param (str) by: The formula of the requested chemical
    or (list) An pair of one cation and one anion for chemtype = 'solid'

        loc(['Ag+', 'Cl-'], 'complex') -> {'AgCl': ..., 'AgCl2-': ..., ...}

    :param (str) chemtype: Leaves None to search all table sequentially. Or use one of the legal type (TABLES)
    to search faster
    :param (bool) clean: Reformat output to python objects.
        E.g: '-3, 1.99' -> {1: -3, 2: 1.99}
        E.g: 'H2SO4, HSO4-, SO42-' -> ['H2SO4', 'HSO4-', 'SO42-']
    :return: (dict) data of the requested chemical. Raises KeyError if nothing was found. [datadict]
        Alternatively, if list searching is selected and chemtype is complex, the query will return multiple matching
        compounds. As such, the datadict will be returned as {[formula1]: [datadict1], [formula2]: [datadict2],...}
    """
    path_db = os.path.join(conf.Paths.PATH_DB_SQL, conf.Paths.FILENAME_DB_SQL)
    tables = conf.IonicEquilibriumOptions.TABLES
    conn = sqlite3.connect(path_db)
    cursor = conn.cursor()

    if isinstance(by, list) and chemtype is None:  # If the user inputted a list, assume chemtype to be solid
        chemtype = 'solid'

    # Checks for illegal type
    if (chemtype is not None) and (chemtype not in tables):
        raise KeyError(f'Illegal type {chemtype}. Either leave chemtype as None to search all tables or \
        it must be a legal type: {tables}')

    if chemtype == 'solid' and isinstance(by, list) and len(by) != 2:
        msg = f'To search for solids/complex using list, \
        please input a list of one cation and one anion.\nYour input: {by}'
        raise TypeError(msg)

    # Search for the chemical if type is known
    compound = []
    if chemtype == 'solid' and isinstance(by, list):  # List search for solid
        query = f'SELECT * FROM solid WHERE cation="{by[0]}" COLLATE NOCASE AND anion="{by[1]}" COLLATE NOCASE '
        cursor.execute(query)
        compound = cursor.fetchall()
        if len(compound) != 1:
            query = f'SELECT * FROM solid WHERE cation="{by[1]}" COLLATE NOCASE AND anion="{by[0]}" COLLATE NOCASE '
            cursor.execute(query)
            compound = cursor.fetchall()
    elif chemtype == 'complex' and isinstance(by, list):  # List search for all complexes that match this criteria
        query = f'SELECT * FROM complex WHERE cation="{by[0]}" COLLATE NOCASE AND ligand="{by[1]}" COLLATE NOCASE '
        cursor.execute(query)
        complexes = cursor.fetchall()
        datadicts = {}
        for cp in complexes:
            datadicts[cp[0]] = dict()
            for i, complex_data in enumerate(cp):
                datadicts[cp[0]][cursor.description[i][0]] = complex_data
            datadicts[cp[0]]['type'] = 'complex'
            datadicts[cp[0]] = clean_data(datadicts[cp[0]])
        return datadicts
    elif chemtype is not None:
        query = f'SELECT * FROM {chemtype} WHERE formula="{by}" COLLATE NOCASE'
        cursor.execute(query)
        compound = cursor.fetchall()
    else:
        # Search for the chemical if type is not known, iterating all legal types until found
        for legal_type in tables:
            query = f'SELECT * FROM {legal_type} WHERE formula="{by}" COLLATE NOCASE'
            cursor.execute(query)
            compound = cursor.fetchall()
            if len(compound) != 0:
                chemtype = legal_type
                break

    # Check if there are more than one compound found or none found. For debugging purposes.
    if not isinstance(by, list) and not isinstance(by, tuple):
        try:
            assert len(compound) == 1
        except AssertionError:
            msg = f'Search for {by} failed.'
            raise KeyError(msg)

    if len(compound) == 0:
        return {}

    # Get the only compound object from the returned searches.
    # Then convert it into a dictionary and add keys (since sqlite SELECT doesn't include names)
    compound = compound[0]
    # Add type to datadict as well
    datadict = {'type': chemtype}
    for i in range(len(compound)):
        datadict[cursor.description[i][0]] = compound[i]
    datadict = clean_data(datadict) if clean else datadict
    conn.commit()
    conn.close()
    return datadict


def iloc(by: str, i: int, chemtype: str = None, clean: bool = True) -> dict:
    """
    Locates a chemical within a database given its formula and type (optional) and return its info in a dictionary.
    With iloc, specify i to identify the conjugate state of the species. (0: conjugate acid, i: conjugate base)
    :param (str) by: The formula of the requested chemical
    :param (str) i: The conjugate state (dissolution tier) of the species. (0: conjugate acid, i: conjugate base)
    :param (str) chemtype: Leaves None to search all table sequentially. Use
    ['anion', 'cation', 'acid', 'base', 'complex'] to search faster.
    :param (bool) clean: Reformat output to python objects.
        E.g: '-3, 1.99' -> {1: -3, 2: 1.99}
        E.g: 'H2SO4, HSO4-, SO42-' -> ['H2SO4', 'HSO4-', 'SO42-']
    :return: (dict) data of the requested chemical. Raises KeyError if nothing was found. [datadict]
    """
    path_db = os.path.join(conf.Paths.PATH_DB_SQL, conf.Paths.FILENAME_DB_SQL)
    tables = conf.IonicEquilibriumOptions.TABLES
    conn = sqlite3.connect(path_db)
    cursor = conn.cursor()

    # Checks for illegal type
    if (chemtype is not None) and (chemtype not in tables):
        raise KeyError(f'Illegal type {chemtype}. Either leave chemtype as None to search all tables or \
        it must be a legal type: {tables}')

    # i can't be anything other can integer
    if not isinstance(i, int):
        msg = f'i must be an integer, not ({type(i)}) {i}'
        raise TypeError(msg)
    # i can't be negative
    if i < 0:
        msg = f'Conjugation state can\'t be a negative value.\nYour input: by={by}, i={i}, chemtype={chemtype}'
        raise ValueError(msg)

    # Search for the chemical if type is known
    compound = []
    if chemtype is not None:
        query = f'SELECT * FROM {chemtype} WHERE formula="{by}" COLLATE NOCASE'
        cursor.execute(query)
        compound = cursor.fetchall()
    else:
        # Search for the chemical if type is not known, iterating all legal types until found
        for legal_type in tables:
            query = f'SELECT * FROM {legal_type} WHERE formula="{by}" COLLATE NOCASE'
            cursor.execute(query)
            compound = cursor.fetchall()
            if len(compound) != 0:
                chemtype = legal_type
                break

    # Check if there are more than one compound found or none found. For debugging purposes.
    try:
        assert len(compound) == 1
    except AssertionError:
        msg = f'Search for {by} failed.'
        raise KeyError(msg)
    # Get the only compound object from the returned searches.
    # Then convert it into a dictionary and add keys (since sqlite SELECT doesn't include names)
    compound = compound[0]
    # Add type to datadict as well
    datadict = {'type': chemtype}
    for n in range(len(compound)):
        datadict[cursor.description[n][0]] = compound[n]
    datadict = clean_data(datadict) if clean else datadict

    # Solid chemtype is illegal
    if chemtype == 'solid':
        msg = 'chemtype "solid" does not have any dissolution tier.'
        raise ValueError(msg)

    # Get the formula of the ith-state conjugated form of the species we're looking for
    try:
        by_i = datadict[conf.IonicEquilibriumOptions.SPECIES_NAME][i]
    except IndexError as err:
        msg = f'Compound {by} does not have a dissolution tier i = {i}'
        raise IndexError(msg)
    conn.commit()
    conn.close()
    return loc(by_i)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='dbmanage.py', description=f'Manage the database for {conf.Paths.APP_NAME}.\n\
    Utilizes 2 files: 1 is an excel file for human editing, 1 is an sql file for computer processing.')

    parser.add_argument('action', type=str,
                        help='dbmake -> convert the excel file to an sql database')

    parser.add_argument('-ip', '--inputpath', type=str,
                        help='The path to the .xlsx file used to create the SQL database.')
    parser.add_argument('-in', '--inputname', type=str, help='The excel file\'s name.')
    parser.add_argument('-op', '--outputpath', type=str, help='The path to the output .xlsx file.')
    parser.add_argument('-on', '--outputname', type=str, help='The SQL file\'s name.')

    args = parser.parse_args()

    if args.action == 'dbmake':
        try:
            input_path = conf.Paths.PATH_DB_XLSX
            input_name = conf.Paths.FILENAME_DB_XLSX
            output_path = conf.Paths.PATH_DB_SQL
            output_name = conf.Paths.FILENAME_DB_SQL
            if args.inputpath:
                input_path = args.inputpath
            if args.inputname:
                input_name = args.inputname
            if args.outputpath:
                output_path = args.outputpath
            if args.outputname:
                output_name = args.outputname
            make_db(input_path, input_name, output_path, output_name)
            print(f'Database creation successful!\nNew file created at {os.path.join(output_path, output_name)}')
        except Exception as e:
            traceback.print_exc()