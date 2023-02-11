import pandas as pd
import sqlite3
import sys
import warnings
from pprint import pprint

import sympy

DB_PATH = 'D:/Coding/Equillum/sentri/db/chemicals.db'


def make_db(update: str = True, xlsx_path: str = 'chemicals.xlsx', db_path: str = DB_PATH):
    """
    Create a .db file with all chemicals information, converted from a .xlsx file.
    :param (bool) update: Set to true to update all chemicals information.
    :param (str) xlsx_path: The path to the xlsx file used to create the .db file.
    :param (str) db_path: Where to place the newly created .db file, shouldn't be changed.
    :return: None
    """
    conn = sqlite3.connect(db_path)
    # Create the three tables cation, anion and acid
    conn.execute(
        '''CREATE TABLE IF NOT EXISTS cation
        (formula           TEXT,
        element            TEXT,
        name               TEXT,
        conjugate     TEXT,
        conjugate_info     TEXT,
        conjugate_count    INTEGER,
        charge             INTEGER,
        multiplier         INTEGER
        );''')
    conn.execute(
        '''CREATE TABLE IF NOT EXISTS complex
        (formula           TEXT,
        cation             TEXT,
        ligand             TEXT,
        info               TEXT,
        charge             INT,
        cation_count       INTEGER,
        ligand_count        INTEGER,
        complex_count      INTEGER
        );''')
    conn.execute(
        '''CREATE TABLE IF NOT EXISTS solid
        (formula           TEXT,
        name               TEXT,
        cation             TEXT,
        anion              TEXT,
        info               TEXT,
        cation_count       INTEGER,
        cation_charge      INTEGER,
        anion_count        INTEGER,
        anion_charge       INTEGER
        );''')
    conn.execute(
        '''CREATE TABLE IF NOT EXISTS anion
        (formula           TEXT,
        element            TEXT,
        name               TEXT,
        conjugate          TEXT,
        conjugate_info     TEXT,
        conjugate_count    INTEGER,
        charge             INTEGER,
        proton_count       INTEGER,
        multiplier         INTEGER,
        special            TEXT,
        special_info       TEXT
        );''')
    conn.execute(
        '''CREATE TABLE IF NOT EXISTS acid
        (formula           TEXT,
        name               TEXT,
        conjugate_info     TEXT,
        conjugate_lookup   TEXT,
        conjugate_count    INTEGER,
        charge             INTEGER,
        multiplier         INTEGER
        );''')

    # Update info if update=True
    if update:
        chemdb_cation = pd.read_excel('chemicals.xlsx', sheet_name='cation')
        chemdb_cation.to_sql('cation', conn, index=False, if_exists='replace')
        chemdb_complex = pd.read_excel('chemicals.xlsx', sheet_name='complex')
        chemdb_complex.to_sql('complex', conn, index=False, if_exists='replace')
        chemdb_solid = pd.read_excel('chemicals.xlsx', sheet_name='solid')
        chemdb_solid.to_sql('solid', conn, index=False, if_exists='replace')
        chemdb_anion = pd.read_excel('chemicals.xlsx', sheet_name='anion')
        chemdb_anion.to_sql('anion', conn, index=False, if_exists='replace')
        chemdb_acid = pd.read_excel('chemicals.xlsx', sheet_name='acid')
        chemdb_acid.to_sql('acid', conn, index=False, if_exists='replace')
    conn.commit()
    conn.close()


def loc(formula: str, chemtype: str = None, db_path: str = DB_PATH) -> dict:
    """
    Locates a chemical within a database given its formula and type (optional) and return its info in a dictionary.
    :param (str) formula: The formula of the requested chemical
    :param (str) chemtype: Leaves None to search all table sequentially. Use ['anion', 'cation', 'acid', 'complex'] to search faster
    :return: (dict) info of the requested chemical. Returns None if nothing was found. [infodict]
    """
    legal_types = ['anion', 'cation', 'acid', 'complex']
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Checks for illegal type
    if (chemtype is not None) and (chemtype not in legal_types):
        raise ValueError(f'Illegal type {chemtype}. Type must be None or in ["anion", "cation", "acid"]')
    # Search for the chemical if type is known
    if chemtype is not None:
        query = f'SELECT * FROM {chemtype} WHERE formula="{formula}"'
        cursor.execute(query)
        compound = cursor.fetchall()
    else:
        # Search for the chemical if type is not known, iterating all legal types until found
        for legal_type in legal_types:
            query = f'SELECT * FROM {legal_type} WHERE formula="{formula}"'
            cursor.execute(query)
            compound = cursor.fetchall()
            if len(compound) != 0:
                chemtype = legal_type
                break
    # Check if there are more than one compound found or none found. For debugging purposes.
    if len(compound) != 1:
        if len(compound) == 0:
            return None
        print('Found', compound)
        raise ValueError(f'Something went wrong while searching for {formula} in the database.')

    # Get the only compound object from the returned searches. Then convert it into a dictionary and add keys (since sqlite SELECT doesn't include names)
    compound = compound[0]
    # Add type to infodict as well
    t = {'type': chemtype}
    for i in range(len(compound)):
        t.update({cursor.description[i][0]:compound[i]})
    conn.commit()
    conn.close()
    return t


def loc_solid(formula: str = None, cation: str = None, anion: str = None, db_path: str = DB_PATH):
    if formula is None and cation is None and anion is None:
        raise ValueError
    if (cation is not None and anion is None) or (cation is not None and anion is None):
        raise ValueError
    if formula is None:
        query = f'SELECT * FROM solid WHERE cation="{cation}" and anion="{anion}"'
    else:
        query = f'SELECT * FROM solid WHERE formula="{formula}"'
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(query)
    compound = cursor.fetchall()
    if len(compound) != 1:
        if len(compound) == 0:
            return None
        print('Found', compound)
        raise ValueError(f'Something went wrong while searching for {formula}({cation}, {anion}) in the database.')
    # Get the only compound object from the returned searches. Then convert it into a dictionary and add keys (since sqlite SELECT doesn't include names)
    compound = compound[0]
    # Add type to infodict as well
    t = {'type': 'solid'}
    for i in range(len(compound)):
        t.update({cursor.description[i][0]: compound[i]})
    conn.commit()
    conn.close()
    return t


def loc_many(complex_pair: list = [], db_path: str = DB_PATH) -> dict:
    """
    Locates many complexes given a list of one cation and one ligand.
    :param (list) complex_pair: A list of exactly one cation (first) and one anion (second)
        loc(chemtype='complex', complex_pair=['Ag', 'Cl-'])
        {'AgCl': ..., 'AgCl2-': ..., ...}
    :param db_path:
    :return: (dict) {[formula1]: [infodict], [formula2]: [infodict],...}
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    assert len(complex_pair) == 2
    # Search for all complexes that match this criteria
    query = f'SELECT * FROM complex WHERE cation="{complex_pair[0]}" AND ligand="{complex_pair[1]}"'
    cursor.execute(query)
    complexes = cursor.fetchall()
    # Convert it into a dictionary and add keys (since sqlite SELECT doesn't include names)
    # Add type to infodict as well
    t = {}
    for complex in complexes:
        t[complex[0]] = dict()
        for i, complex_info in enumerate(complex):
            t[complex[0]][cursor.description[i][0]] = complex_info
        t[complex[0]]['type'] = 'complex'
    conn.commit()
    conn.close()
    return t


def iloc(formula: str, i: int, chemtype: str = None, db_path: str = DB_PATH) -> dict:
    """
    Locates a chemical within a database given its formula and type (optional) and return its info in a dictionary.
    With iloc, specify i to identify the conjugate state of the species. (0: extreme acidic, n: extreme basic)
    :param (str) formula: The formula of the requested chemical
    :param (str) type: Leaves None to search all table sequentially. Use ['anion', 'cation', 'acid'] to search faster
    :return: (dict) info of the requested chemical. Returns None if nothing was found.
    """
    legal_types = ['anion', 'cation', 'acid']
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Checks for illegal type
    if (chemtype is not None) and (chemtype not in legal_types):
        raise ValueError(f'Illegal type {chemtype}. Type must be None or in ["anion", "cation", "acid"]')

    # Search for the chemical if type is known
    if chemtype is not None:
        query = f'SELECT * FROM {chemtype} WHERE formula="{formula}"'
        cursor.execute(query)
        compound = cursor.fetchall()
    else:
        # Search for the chemical if type is not known, iterating all legal types until found
        for legal_type in legal_types:

            query = f'SELECT * FROM {legal_type} WHERE formula="{formula}"'
            cursor.execute(query)
            compound = cursor.fetchall()
            if len(compound) != 0:
                chemtype = legal_type
                break
    # Check if there are more than one compound found or none found. For debugging purposes.
    if len(compound) != 1:
        if len(compound) == 0:
            warnings.warn(f'Unable to locate {formula} in the database.')
            return None
        print('Found', compound)
        raise ValueError(f'Something went wrong while searching for {formula} in the database.')

    # Get the only compound object from the returned searches. Then convert it into a dictionary and add keys (since sqlite SELECT doesn't include names)
    compound = compound[0]
    # Add type to infodict as well
    t = {'type': chemtype}
    for n in range(len(compound)):
        t.update({cursor.description[n][0]:compound[n]})
    # Get the formula of the ith-state conjugated form of the species we're looking for
    formula = t['conjugate'].splitlines()[i]
    conn.commit()
    conn.close()
    return loc(formula)


def loc_elem(element: str, db_path: str = DB_PATH) -> dict:
    """
    Locates an element within a database given its formula and type (optional) and return its info in a dictionary.
    :param (str) formula: The formula of the requested chemical
    :param (str) chemtype: Leaves None to search all table sequentially. Use ['anion', 'cation', 'acid', 'complex'] to search faster
    :return: (dict) info of the requested chemical. Returns None if nothing was found. [infodict]
    """
    legal_types = ['anion', 'cation']
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Search for the chemical if type is not known, iterating all legal types until found
    for legal_type in legal_types:
        query = f'SELECT * FROM {legal_type} WHERE element="{element}"'
        cursor.execute(query)
        compound = cursor.fetchall()
        if len(compound) != 0:
            chemtype = legal_type
            break
    # Check if there are more than one compound found or none found. For debugging purposes.

    # Get the only compound object from the returned searches. Then convert it into a dictionary and add keys (since sqlite SELECT doesn't include names)
    return_dict = dict()
    for comp in compound:
        t = dict()
        for i in range(len(comp)):
            t.update({cursor.description[i][0]:comp[i]})
        return_dict[t['formula']] = t['charge']
    conn.commit()
    conn.close()
    return return_dict


def extract_const(infostr: str) -> dict:
    """
    Given an infostring on conjugate (e.g: H2AsO3-, pKa1 = 9.29), this function extracts all constants and return a dict
    :param (str) infostr: The string containing information
    :return: (dict)
    """
    const = {}
    # Return empty list if it is a special type (infostr is None)
    if infostr is None:
        return dict()
    # Seperate at each line for each conjugate forms
    for conjugate_tier in infostr.splitlines():
        current_tier_const = float(conjugate_tier.split(',')[-1].split('=')[-1].strip())
        const.update({len(const) + 1 : current_tier_const})
    return const


def extract_conjugate(infostr: str) -> dict:
    """
    Given an infostring on conjugate (e.g: H2AsO3-, pKa1 = 9.29), this function extracts all constants and return a dict
    :param (str) infostr: The string containing information
    :return: (dict)
    """
    conj = {}
    # Return empty list if it is a special type (infostr is None)
    if infostr is None:
        return dict()
    # Seperate at each line for each conjugate forms
    for conjugate_tier in infostr.splitlines():
        current_tier_formula = conjugate_tier.split(',')[0].strip()
        current_tier_const = float(conjugate_tier.split(',')[-1].split('=')[-1].strip())
        conj.update({current_tier_formula: current_tier_const})
    return conj


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


def parse(origin: str, to: str, obj):
    """
    Convert an object from one type to another. Check documentation for allowed conversion
        1) 'sympy.Expr' -> 'wolfram'     2) 'sympy.Eq' -> wolfram'
        3) 'sympy.Expr' -> 'matlab'      4) 'sympy.Eq' -> 'matlab
    :param obj: The object that needs converting
    :return:
    """
    return_obj = ''
    def sympyexpr_to_wolfram(obj):
        obj = str(obj)
        obj = obj.replace("**", "^").replace("e-", "*^-")
        # sqrt conversion
        for i, char in enumerate(obj):
            if char == 's' and obj[i:i + 5] == 'sqrt(':
                # Entering replacement mode
                open_bracket = 1
                for n, next_char in enumerate(obj[i + 5:]):
                    n = n + i + 5
                    if next_char == '(':
                        open_bracket += 1
                    if next_char == ')':
                        open_bracket -= 1
                    if open_bracket == 0:
                        obj = obj[:n] + ']' + obj[n + 1:]
                        break
                obj = obj[:i] + 'Sqrt[' + obj[i + 5:]
        return obj
    if origin == 'sympy.Eq' and not isinstance(obj, sympy.Eq):
        raise Exception
    if origin == 'sympy.Expr' and to == 'wolfram':
        return sympyexpr_to_wolfram(obj)
    if origin == 'sympy.Eq' and to == 'wolfram':
        eq_l = sympyexpr_to_wolfram(obj.lhs)
        eq_r = sympyexpr_to_wolfram(obj.rhs)
        return f'{eq_l} == {eq_r}'
    if origin == 'sympy.Expr' and to == 'matlab':
        return str(obj).replace('**','.^')
    if origin == 'sympy.Eq' and to == 'matlab':
        return f"{str(obj.lhs).replace('**','.^')} == {str(obj.rhs).replace('**','.^')}"
    raise Exception


def inspect_formula(formula: str) -> dict:
    cation = dict()
    anion = dict()
    neutral = dict()
    return_dict = {'cation': cation, 'anion': anion, 'neutral': neutral}
    elements = dict()
    indices = dict()
    # Creates a dictionary of all elements and their respective count in the formula (without multiplier)
    # Fe2(SO4)3 -> {'Fe': 2, 'S': 1, 'O': 4}
    for i, char in enumerate(formula):
        if char.isupper():
            curr_elem = char
            curr_elem_count = 0
            curr_elem_index = i
            for next_char in formula[i+1:]:
                # If the next character after an element is another element (upper alpha), set count to 1 and end
                # If next char is still that element (lower alpha), add it to the formula string and continue
                # If next char is its count explicitly declared (more than 1), set it to that number and continue
                # (continue above) If another number follow it, use this formula to get count: curr_elem_count * 10 + int(next_char)
                if next_char.isalnum():
                    if next_char.islower():
                        curr_elem += next_char
                    elif next_char.isupper():
                        curr_elem_count = 1
                        break
                    elif next_char.isnumeric():
                        if curr_elem_count != 0:
                            curr_elem_count = curr_elem_count * 10 + int(next_char)
                            continue
                        curr_elem_count = int(next_char)
                    else:
                        break
                else:
                    break

            elements[curr_elem] = curr_elem_count
            indices[curr_elem_index] = curr_elem
    # Creates a dictionary of all elements and their multiplier (if any) defined by the number after the brackets.
    multiplier = {key: 1 for key in elements.keys()}
    curr_multiplier = 1
    in_round_bracket = False
    round_bracket_count = 0
    round_bracket_start = 0
    round_bracket_end = 0
    for i, char in enumerate(formula):
        if char == '(':
            in_round_bracket = True
            round_bracket_start = i
            round_bracket_count += 1
        if char == ')':
            round_bracket_count -= 1
            assert round_bracket_count >= 0
        if in_round_bracket and round_bracket_count == 0:
            round_bracket_end = i
            curr_multiplier = int(formula[i+1])
            in_round_bracket = False
            for n in range(round_bracket_start, round_bracket_end):
                if indices.get(n) is not None:
                    multiplier[indices[n]] = curr_multiplier
    elements = {key: val * multiplier[key] for key, val in elements.items()}

    print(elements)
    print(multiplier)
    print('')

    print(loc_elem([elem for elem in elements.keys()][0]))


if __name__ == '__main__':
    args = set(sys.argv)
    if ('make_db' in args) or ('makedb' in args) or ('create_db' in args) or ('dbmake' in args):
        try:
            make_db()
        except Exception as e:
            print(e)
        print('Database creation sucessful!\nNew file created at sentri/db/chemicals.db')
