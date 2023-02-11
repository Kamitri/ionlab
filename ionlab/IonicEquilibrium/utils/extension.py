import sympy
from gekko import GEKKO


class SympyExtension:

    @staticmethod
    def to_python(obj: sympy.Eq, name: str):
        return 3

    @staticmethod
    def to_gekko_model(eqns: list[sympy.Eq] | dict[str, sympy.Eq] | sympy.Eq,
                       variables: dict[str, float],
                       constants: dict[str, float]) -> list:
        """
        Convert a list of sympy equations into a GEKKO model ready for solving
        :param eqns: A list of sympy equations / one sympy equation / a dictionary of sympy equations
        :param variables: A dictionary with the ion symbol (name) as key and its estimation as value
            Example: {'h': 10 ** -14, 'xani0st0': 10 ** -14}
        :param constants: A dictionary with the constant symbol (name) as key and its numeric value as value
            Example: {'pKaani0st1': -3, 'pKaani0st2': 1.99}
        :return: A list consisting 2 objects. The former is the GEKKO model. The latter is the variable for accessing
        the solutions.
        """
        if isinstance(eqns, sympy.Eq):
            eqns = [eqns]
        if isinstance(eqns, dict):
            eqns = [eqn for eqn_name, eqn in eqns.items()]

        m = GEKKO()
        constants = {const_name: m.Const(float(const_value)) for const_name, const_value in constants.items()}
        variables = {var_name: m.Var(float(var_value)) for var_name, var_value in variables.items()}
        # Set lower bounds of all variables to be at least 0
        for variable in variables.values():
            variable.LOWER = 0

        gekko_eqns = list()
        for eqn in eqns:
            eqn_str = str(eqn.lhs) + ' == ' + str(eqn.rhs)
            gekko_eqns.append(eval(eqn_str, {"__builtins__": None}, {**constants, **variables, 'eqn': eqn}))
        m.Equations(gekko_eqns)
        return [m, variables]

    @staticmethod
    def ensure_pos(expr: sympy.Expr, assumption: int = 1) -> bool:
        """
        Given a sympy expression, this function will return True if this expression evaluates to larger than 0 and otherwise.
        :param (sympy.Expr) expr: Sympy expression to be evaluated.
        :param (int) assumption: This assumes that all free variables within the expression will be equal to this value.
        :return: (bool)
        """
        free_symbols = expr.free_symbols
        subs_dict = {symbol: 1 for symbol in free_symbols}
        evaluation = expr.evalf(subs=subs_dict)
        if evaluation > 0:
            return True
        return False


if __name__ == '__main__':
    # m = GEKKO()
    # h = m.Var(value=10e-14)
    # xani0st0 = m.Var(value=10e-14)
    # pKaani0st1 = m.Const(-3)
    # pKaani0st2 = m.Const(1.99)
    h = sympy.Symbol('h')
    xani0st0 = sympy.Symbol('xani0st0')
    pKaani0st1 = sympy.Symbol('pKaani0st1')
    pKaani0st2 = sympy.Symbol('pKaani0st2')
    obj1 = sympy.Eq(h, 1.0e-14 / h + xani0st0 / (10 ** pKaani0st1 * h) + 2 * xani0st0 / (
                10 ** pKaani0st1 * 10 ** pKaani0st2 * h ** 2))
    obj2 = sympy.Eq(0.06, xani0st0 + xani0st0 / (10 ** pKaani0st1 * h) + xani0st0 / (
                10 ** pKaani0st1 * 10 ** pKaani0st2 * h ** 2))
    consts = {'pKaani0st1': -3, 'pKaani0st2': 1.99}
    eq_vars = {'h': 10 ** -14, 'xani0st0': 10 ** -14}
    model, model_vars = SympyExtension.to_gekko_model(eqns=[obj1, obj2], variables=eq_vars, constants=consts)
    model.solve(disp=False)
    print(model_vars['h'].value[0], model_vars)
