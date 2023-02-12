from gekko import GEKKO
import sympy
from typing import Union, List, Dict
from ionlab.IonicEquilibrium.utils.extension import SympyExtension as ext
import numpy as np


class GEKKOSolver:
    def __init__(self, strict_check: bool = True):
        # TODO: m.options.SOLVERS allow for more algorithms
        self.algorithms = {
            'Algorithm 1-STRICT.APPROX_S': self.algorithm1,
            'Algorithm 2-STRICT.APPROX_S': self.algorithm2,
            'Algorithm 3-STRICT.APPROX_S': self.algorithm3,
            'Algorithm 1.APPROX_S': self.algorithm1,
            'Algorithm 1-STRICT.APPROX_M': self.algorithm1,
            'Algorithm 2-STRICT.APPROX_M': self.algorithm2,
            'Algorithm 3-STRICT.APPROX_M': self.algorithm3,
            'Algorithm 1.APPROX_M': self.algorithm1,
            'Algorithm 1-STRICT.APPROX_L': self.algorithm1,
            'Algorithm 2-STRICT.APPROX_L': self.algorithm2,
            'Algorithm 3-STRICT.APPROX_L': self.algorithm3,
            'Algorithm 1.APPROX_L': self.algorithm1,
            'Algorithm 2.APPROX_S': self.algorithm2,
            'Algorithm 3.APPROX_S': self.algorithm3,
            'Algorithm 2.APPROX_M': self.algorithm2,
            'Algorithm 3.APPROX_M': self.algorithm3,
            'Algorithm 2.APPROX_L': self.algorithm2,
            'Algorithm 3.APPROX_L': self.algorithm3,
        }
        self.eq_molarity_dict = {}
        self.consts = {}
        self.strict_check = strict_check
        self.approximations = {'APPROX_S': 10**-14, 'APPROX_M': 10**-7, 'APPROX_L': 0.1}

    def compute_eq(self, eqns: List[sympy.Eq], variables: List[str], constants: Dict[str, float]) -> dict:
        self.consts = constants
        try:
            next_algorithm_key = next(iter(self.algorithms))
        except Exception as e:
            print('No algorithms were able to compute the equilibrium for this system.')
            print(e)
            return {}
        try:
            strict_mode = True if 'STRICT' in next_algorithm_key else False
            index = next_algorithm_key.find('APPROX_')
            approximation_lv = self.approximations[next_algorithm_key[index:]]
            self.algorithms[next_algorithm_key]\
                (eqns=eqns, variables=variables, constants=constants,
                 strict_mode=strict_mode, approximation=approximation_lv)
            if not self.output_check(eqns, self.eq_molarity_dict):
                raise Exception
        except Exception as e:
            # print(f'{next_algorithm_key} failed.')
            # print(e)
            self.algorithms.pop(next_algorithm_key)
            return self.compute_eq(eqns, variables, constants)
        return self.eq_molarity_dict

    def algorithm1(self, eqns: List[sympy.Eq], variables: List[str],
                   constants: Dict[str, float], approximation: float, strict_mode=True) -> None:
        variables = {var: approximation for var in variables}
        model, variables = ext.to_gekko_model(eqns=eqns, variables=variables, constants=constants)
        if strict_mode:
            model.options.RTOL = 10e-10
        model.options.SOLVER = 1
        model.solve(disp=False)
        self.eq_molarity_dict = {key: val.value[0] for key, val in variables.items()}
        self.eq_molarity_dict['pH'] = -np.log10(self.eq_molarity_dict['h'])

    def algorithm2(self, eqns: List[sympy.Eq], variables: List[str],
                   constants: Dict[str, float], approximation: float, strict_mode=True) -> None:
        variables = {var: approximation for var in variables}
        model, variables = ext.to_gekko_model(eqns=eqns, variables=variables, constants=constants)
        if strict_mode:
            model.options.RTOL = 10e-10
        model.options.SOLVER = 2
        model.solve(disp=False)
        self.eq_molarity_dict = {key: val.value[0] for key, val in variables.items()}
        self.eq_molarity_dict['pH'] = -np.log10(self.eq_molarity_dict['h'])

    def algorithm3(self, eqns: List[sympy.Eq], variables: List[str],
                   constants: Dict[str, float], approximation: float, strict_mode=True) -> None:
        variables = {var: approximation for var in variables}
        model, variables = ext.to_gekko_model(eqns=eqns, variables=variables, constants=constants)
        if strict_mode:
            model.options.RTOL = 10e-10
        model.options.SOLVER = 3
        model.solve(disp=False)
        self.eq_molarity_dict = {key: val.value[0] for key, val in variables.items()}
        self.eq_molarity_dict['pH'] = -np.log10(self.eq_molarity_dict['h'])

    def output_check(self, eqns: List[sympy.Eq], eq_molarity_dict: dict) -> bool:
        # Check if any output is negative
        if not all([val > 0 for val in eq_molarity_dict.values()]):
            # print('Computation failed.')
            return False
        tol = 10e-8 if self.strict_check else 10e-4
        for eqn in eqns:
            lhs_val = eqn.lhs.subs({**eq_molarity_dict, **self.consts})
            rhs_val = eqn.rhs.subs({**eq_molarity_dict, **self.consts})
            if abs(sympy.N(lhs_val - rhs_val, 15)) > tol:
                return False
        return True