from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wl, wlexpr
from wolframclient.language.expression import WLFunction, WLSymbol
import ionlab.conf as conf
import sympy
from ionlab.IonicEquilibrium.db import process


class WolframSolver:
    def __init__(self, debug: bool = False):
        self.algorithms = {
            'Algorithm 1': self.algorithm1,  # Do not put . in the name as the program use the text after . as args
            'Algorithm 1-RE': self.algorithm1_reroot,
            'Algorithm 2': self.algorithm2,  # Do not put . in the name as the program use the text after . as args
            'Algorithm 2-RE': self.algorithm2_reroot,
            'Algorithm 1-WP.1': self.algorithm1_wp,
            'Algorithm 1-WP.20': self.algorithm1_wp,
            'Algorithm 1-WP.50': self.algorithm1_wp,
            'Algorithm 1-WP.100': self.algorithm1_wp,
            'Algorithm 1-WP.MachinePrecision': self.algorithm1_wp,
            'Algorithm 1-RE-WP.1': self.algorithm1_reroot_wp,
            'Algorithm 1-RE-WP.20': self.algorithm1_reroot_wp,
            'Algorithm 1-RE-WP.50': self.algorithm1_reroot_wp,
            'Algorithm 1-RE-WP.100': self.algorithm1_reroot_wp,
            'Algorithm 1-RE-WP.MachinePrecision': self.algorithm2_reroot_wp,
            'Algorithm 2-RE-WP.1': self.algorithm2_reroot_wp,
            'Algorithm 2-RE-WP.20': self.algorithm2_reroot_wp,
            'Algorithm 2-RE-WP.50': self.algorithm2_reroot_wp,
            'Algorithm 2-RE-WP.100': self.algorithm2_reroot_wp,
            'Algorithm 2-RE-WP.MachinePrecision': self.algorithm2_reroot_wp}
        if conf.Paths.WOLFRAM_KERNEL is None:
            raise ValueError('Wolfram Kernel path hasn\'t been set. Please point conf.Paths.WOLFRAM_KERNEL to your '
                             'WolframKernel.exe')
        self.kernel_path = conf.Paths.WOLFRAM_KERNEL
        self.eq_count = 0
        self.syms_to_solve = list()
        self.debug = debug
        with open('wolfram_solver.m', 'w') as wolfram_solver:
            wolfram_solver.write('')
        self.session = WolframLanguageSession(self.kernel_path)

    def log(self, text: str):
        with open('wolfram_solver.m', 'a') as wolfram_solver:
            if text[-2:] != '\n':
                text += '\n'
            wolfram_solver.write(text)

    def write_expression(self, symbol: str, value: str):
        expr = f'{symbol} = {value};'.replace('e-', '*^-')
        self.session.evaluate(expr)
        if self.debug:
            self.log(expr)

    def write_equation(self, eq: sympy.Eq):
        # Parse sympy.Eq into a valid wolfram equation (the Sqrt[], *^ things...)
        wolfram_eq = f'eq{self.eq_count} = {process.parse(origin="sympy.Eq", to="wolfram", obj=eq)}'
        self.session.evaluate(wolfram_eq)
        if self.debug:
            self.log(wolfram_eq)
        self.eq_count += 1

    def parse_output(self, wlrule: tuple) -> dict:
        """
        Converts a Mathematica rule to a python dictionary
        E.g: (Rule[Global`h, 0.1], Rule[Global`xcat0st0, 0.01]) -> {'h': 0.1, 'xcat0st0': 0.01}
        :param (tuple) wlrule: Mathematica rule
        :return:
        """
        pydict = {}
        for symbol, value in wlrule:
            pydict[str(symbol).split('`')[-1]] = value
        return pydict

    def compute_eq(self) -> dict:
        try:
            next_algorithm_key = next(iter(self.algorithms))
        except Exception as e:
            print('No algorithms were able to compute the equilibrium for this system.')
            print(e)
            return {}
        try:
            if '.' in next_algorithm_key:
                working_precision = next_algorithm_key.split('.')[-1]
                self.algorithms[next_algorithm_key](working_precision)
            self.algorithms[next_algorithm_key]()
            eq_molarity_dict = self.session.evaluate('answer')
            eq_molarity_dict = self.parse_output(eq_molarity_dict)
            if not self.output_check(eq_molarity_dict):
                raise Exception
        except Exception as e:
            # print(f'{next_algorithm_key} failed.')
            self.algorithms.pop(next_algorithm_key)
            return self.compute_eq()
        self.session.terminate()
        return eq_molarity_dict

    @staticmethod
    def output_check(eq_molarity_dict) -> bool:
        """
        Checks the eq_molarity dict if all values are real number and greater than 0
        :param eq_molarity_dict:
        :return:
        """
        for ion_symbol, molarity in eq_molarity_dict.items():
            if isinstance(molarity, WLFunction) or isinstance(molarity, WLSymbol):
                return False
            if molarity <= 0:
                return False
        return True

    def algorithm1(self) -> None:
        eqs = ''
        for i in range(self.eq_count):
            eqs = process.join(eqs, f'eq{i}', '&&')
        if self.debug:
            self.log(f'answer = Solve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals][[1]];')
            self.log('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]]')
            self.log('answer')
        self.session.evaluate(
            f'answer = Solve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals][[1]];')
        self.session.evaluate('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]];')

    def algorithm1_reroot(self) -> None:
        eqs = ''
        for i in range(self.eq_count):
            eqs = process.join(eqs, f'eq{i}', '&&')
        syms_reroot = ''
        for i, sym in enumerate(self.syms_to_solve):
            i += 1
            syms_reroot = syms_reroot + '{' + sym + ', ' + f'answer[[{i}]][[2]]' + '}'
            syms_reroot += ', ' if i != len(self.syms_to_solve) else ''
        if self.debug:
            self.log(f'answer = Solve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals][[1]];')
            self.log(f'answer = FindRoot[{eqs}, {{{syms_reroot}}}];')
            self.log('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]]')
            self.log('answer')
        self.session.evaluate(
            f'answer = Solve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals][[1]];')
        self.session.evaluate(f'answer = FindRoot[{eqs}, {{{syms_reroot}}}];')
        self.session.evaluate('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]];')

    def algorithm1_wp(self, working_precision) -> None:
        eqs = ''
        for i in range(self.eq_count):
            eqs = process.join(eqs, f'eq{i}', '&&')
        if self.debug:
            self.log(f'answer = Solve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals, WorkingPrecision->{working_precision}][[1]];')
            self.log('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]]')
            self.log('answer')
        self.session.evaluate(
            f'answer = Solve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals, WorkingPrecision->{working_precision}][[1]];')
        self.session.evaluate('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]];')

    def algorithm1_reroot_wp(self, working_precision) -> None:
        eqs = ''
        for i in range(self.eq_count):
            eqs = process.join(eqs, f'eq{i}', '&&')
        syms_reroot = ''
        for i, sym in enumerate(self.syms_to_solve):
            i += 1
            syms_reroot = syms_reroot + '{' + sym + ', ' + f'answer[[{i}]][[2]]' + '}'
            syms_reroot += ', ' if i != len(self.syms_to_solve) else ''
        if self.debug:
            self.log(f'answer = Solve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals, WorkingPrecision->{working_precision}][[1]];')
            self.log(f'answer = FindRoot[{eqs}, {{{syms_reroot}}}];')
            self.log('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]]')
            self.log('answer')
        self.session.evaluate(
            f'answer = Solve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals, WorkingPrecision->{working_precision}][[1]];')
        self.session.evaluate(f'answer = FindRoot[{eqs}, {{{syms_reroot}}}];')
        self.session.evaluate('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]];')

    def algorithm2(self) -> None:
        eqs = ''
        for i in range(self.eq_count):
            eqs = process.join(eqs, f'eq{i}', '&&')
        if self.debug:
            self.log(
                f'answer = NSolve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals][[1]];')
            self.log('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]]')
            self.log('answer')
        self.session.evaluate(
            f'answer = NSolve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals][[1]];')
        self.session.evaluate('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]];')

    def algorithm2_reroot(self) -> None:
        eqs = ''
        for i in range(self.eq_count):
            eqs = process.join(eqs, f'eq{i}', '&&')
        syms_reroot = ''
        for i, sym in enumerate(self.syms_to_solve):
            i += 1
            syms_reroot = syms_reroot + '{' + sym + ', ' + f'answer[[{i}]][[2]]' + '}'
            syms_reroot += ', ' if i != len(self.syms_to_solve) else ''
        if self.debug:
            self.log(f'answer = NSolve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals][[1]];')
            self.log(f'answer = FindRoot[{eqs}, {{{syms_reroot}}}];')
            self.log('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]]')
            self.log('answer')
        self.session.evaluate(
            f'answer = NSolve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals][[1]];')
        self.session.evaluate(f'answer = FindRoot[{eqs}, {{{syms_reroot}}}];')
        self.session.evaluate('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]];')

    def algorithm2_wp(self, working_precision) -> None:
        eqs = ''
        for i in range(self.eq_count):
            eqs = process.join(eqs, f'eq{i}', '&&')
        if self.debug:
            self.log(
                f'answer = NSolve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals, WorkingPrecision->{working_precision}][[1]];')
            self.log('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]]')
            self.log('answer')
        self.session.evaluate(
            f'answer = NSolve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals, WorkingPrecision->{working_precision}][[1]];')
        self.session.evaluate('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]];')

    def algorithm2_reroot_wp(self, working_precision) -> None:
        eqs = ''
        for i in range(self.eq_count):
            eqs = process.join(eqs, f'eq{i}', '&&')
        syms_reroot = ''
        for i, sym in enumerate(self.syms_to_solve):
            i += 1
            syms_reroot = syms_reroot + '{' + sym + ', ' + f'answer[[{i}]][[2]]' + '}'
            syms_reroot += ', ' if i != len(self.syms_to_solve) else ''
        if self.debug:
            self.log(f'answer = NSolve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals, WorkingPrecision->{working_precision}][[1]];')
            self.log(f'answer = FindRoot[{eqs}, {{{syms_reroot}}}];')
            self.log('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]]')
            self.log('answer')
        self.session.evaluate(
            f'answer = NSolve[{eqs}, {{{",".join(self.syms_to_solve)}}}, NonNegativeReals, WorkingPrecision->{working_precision}][[1]];')
        self.session.evaluate(f'answer = FindRoot[{eqs}, {{{syms_reroot}}}];')
        self.session.evaluate('AppendTo[answer,pH->-Log10[answer[[1]][[2]]]];')