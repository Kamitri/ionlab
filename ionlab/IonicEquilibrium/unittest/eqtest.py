import unittest
from IonicEquilibrium.db.common_chemicals import *
from IonicEquilibrium.db.matter import Compound, Solution
from IonicEquilibrium.analyze_solution import Lab

HCl = Compound(cation='H+', anion='Cl-', cation_count=1, anion_count=1)
HBr = Compound(cation='H+', anion='Br-', cation_count=1, anion_count=1)
H2SO4 = Compound(cation='H+', anion='SO42-', cation_count=2, anion_count=1)
HCN = Compound(cation='H+', anion='CN-', cation_count=1, anion_count=1)
HClO4 = Compound(formula='HClO4', cation_count=1, anion_count=1)
NaCl = Compound(cation='Na+', anion='Cl-', cation_count=1, anion_count=1)
NaOH = Compound(cation='Na+', anion='OH-', cation_count=1, anion_count=1)
NaCN = Compound(cation='Na+', anion='CN-', cation_count=1, anion_count=1)
NaHS = Compound(cation='Na+', anion='HSO3-', cation_count=1, anion_count=1)
Na3PO4 = Compound(cation='Na+', anion='PO43-', cation_count=3, anion_count=1)
Na2HPO4 = Compound(cation='Na+', anion='HPO42-', cation_count=2, anion_count=1)
NaH2PO4 = Compound(cation='Na+', anion='H2PO4-', cation_count=1, anion_count=1)
BaCl2 = Compound(cation='Ba2+', anion='Cl-', cation_count=1, anion_count=2)
Na2CO3 = Compound(cation='Na+', anion='CO32-', cation_count=2, anion_count=1)
HCOOH = Compound(cation='H+', anion='HCOO-', cation_count=1, anion_count=1)
CH3COOH = Compound(cation='H+', anion='CH3COO-', cation_count=1, anion_count=1)
NH2OH = Compound(formula='NH2OH')
CH3NH3Cl = Compound(cation='CH3NH3+', anion='Cl-', cation_count=1, anion_count=1)
NH3 = Compound(formula='NH3')
NH4Cl = Compound(cation='NH4+', anion='Cl-', cation_count=1, anion_count=1)
KBr = Compound('K+', 'Br-', 1, 1)
KOH = Compound(cation='K+', anion='OH-', cation_count=1, anion_count=1)
Cr2SO43 = Compound(cation='Cr3+', anion='SO42-', cation_count=2, anion_count=3)
K2Cr2O7 = Compound(cation='K+', anion='Cr2O72-', cation_count=2, anion_count=1)
CdBr2 = Compound(cation='Cd2+', anion='Br-', cation_count=1, anion_count=2)
select_solver = 'gekko'


class StrongAcid(unittest.TestCase):
    def test_HCl_01M(self):
        solA = Solution(1)
        solA.add(HCl, init_M=0.1)
        ana = Lab(solA)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        self.assertAlmostEqual(pH, 1, 3)

    def test_HCl_1M(self):
        solA = Solution(1)
        solA.add(HCl, init_M=1)
        ana = Lab(solA)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        self.assertAlmostEqual(pH, 0, 2)

    def test_HCl_00015M(self):
        solA = Solution(1)
        solA.add(HCl, init_M=0.00015)
        ana = Lab(solA)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        self.assertAlmostEqual(pH, 3.824, 3)

    def test_H2SO4_01M(self):
        solA = Solution(1)
        solA.add(H2SO4, init_M=0.1)
        ana = Lab(solA)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        self.assertAlmostEqual(pH, 0.964169, 3)


class WeakAcid(unittest.TestCase):
    pass


class Neutral(unittest.TestCase):
    def test_NaCl_01M(self):
        solA = Solution(1)
        solA.add(NaCl, init_M=0.1)
        ana = Lab(solA)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        self.assertAlmostEqual(pH, 7, 2)

    def test_NaCl_1M(self):
        solA = Solution(1)
        solA.add(NaCl, init_M=1)
        ana = Lab(solA)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        self.assertAlmostEqual(pH, 7, 2)


class NTD_Workbook_Chapter_II(unittest.TestCase):
    def test_II_3_2(self):
        sol1 = Solution(0.001)
        sol1.add(NaOH, init_M=0.2)
        sol2 = Solution(0.001)
        sol2.add(HCl, init_M=0.05)
        sol2.add(CH3COOH, init_M=0.18)
        sol1 = sol2 + sol1
        ana = Lab(sol1)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        self.assertAlmostEqual(pH, 5.45842, 2)

    def test_II_3_3(self):
        sol1 = Solution(1)
        sol1.add(KBr, init_M=0.06)
        sol1.add(K2Cr2O7, init_M=0.01)
        sol1.add(H2SO4, init_M=1)
        sol1.add(Cr2SO43, init_M=0.001)
        ana = Lab(sol1)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        # TODO: Difference here check later
        self.assertAlmostEqual(pH, 0.0016533, 2)

    def test_II_3_6(self):
        sol1 = Solution(0.03 * 10 ** -3)
        sol1.add(NaOH, init_M=0.001)
        sol2 = Solution(100 * 10 ** -3)
        sol2.add(NaCl, init_M=0.1)
        sol1 = sol2 + sol1
        ana = Lab(sol1)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        self.assertAlmostEqual(pH, 7.52, 2)

    def test_II_3_35(self):
        sol1 = Solution(1)
        sol1.add(CH3COOH, init_M=0.01)
        print(sol1.__dict__)
        ana = Lab(sol1)
        print(ana.__dict__)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        self.assertAlmostEqual(pH, 3.39, 2)

    def test_II_3_36(self):
        sol1 = Solution(20*10**-3)
        sol1.add(NH3, init_M=1.5*10**-3)
        sol2 = Solution(40*10**-3)
        sol2.add(HCl, init_M=7.5*10**-4)
        sol1 = sol2 + sol1
        ana = Lab(sol1)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        self.assertAlmostEqual(pH, 6.26, 2)

    def test_II_3_37a(self):
        sol1 = Solution(1)
        sol1.add(NH2OH, init_M=10**-3)
        ana = Lab(sol1)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        # TODO: There's a slight difference here? 8.49
        self.assertAlmostEqual(pH, 8.5, 1)

    def test_II_3_37b(self):
        sol1 = Solution(1)
        sol1.add(NH2OH, init_M=10**-5)
        ana = Lab(sol1)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        # TODO: There's a slight difference here? 7.49
        self.assertAlmostEqual(pH, 7.5, 1)

    def test_II_3_72(self):
        sol1 = Solution(0.4)
        sol1.add(HCl, init_M=2*10**-4)
        sol1.add(NH4Cl, init_M=10 ** -2)
        ana = Lab(sol1)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        self.assertAlmostEqual(pH, 3.70, 2)

    def test_II_3_73(self):
        sol1 = Solution(0.4)
        sol1.add(NaOH, init_M=2.51*10**-2)
        sol1.add(NH4Cl, init_M=2.5*10**-2)
        ana = Lab(sol1)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        self.assertAlmostEqual(pH, 10.85, 2)

    def test_II_3_74(self):
        sol1 = Solution(0.03*10**-3)
        sol1.add(KOH, init_M=0.084)
        sol2 = Solution(0.1)
        sol2.add(HCOOH, init_M=2.45 * 10 ** -5)
        sol1 = sol2 + sol1
        ana = Lab(sol1)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        self.assertAlmostEqual(pH, 7.85, 2)

    def test_II_3_75(self):
        sol1 = Solution(50*10**-3)
        sol1.add(NH3, init_M=2*10**-3)
        sol2 = Solution(50*10**-3)
        sol2.add(H2SO4, init_M=2*10**-3)
        sol1 = sol2 + sol1
        ana = Lab(sol1)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        self.assertAlmostEqual(pH, 3.04, 2)

    def test_II_3_76(self):
        sol1 = Solution(1)
        sol1.add(HCN, init_M=10**-4)
        sol1.add(CH3NH3Cl, init_M=2*10**-3)
        ana = Lab(sol1)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        self.assertAlmostEqual(pH, 6.49, 2)


class NTD_Workbook_Chapter_III(unittest.TestCase):
    def test_III_2_7(self):
        sol1 = Solution(1)
        sol1.add(CdBr2, init_M=0.01)
        sol1.add(HBr, init_M=1)
        ana = Lab(sol1)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        symbol_of_CdBr = ana.__dict__['logbook']['CdBr+']
        self.assertAlmostEqual(ana.__dict__['complex_reference'][symbol_of_CdBr]['value'], 0.0006388, 2)


class NTD_Workbook_Chapter_V(unittest.TestCase):
    def test_V_1_17(self):
        sol1 = Solution(1)
        sol1.add(BaCl2, init_M=1)
        sol1.add(Na2CO3, init_M=1)
        ana = Lab(sol1)
        ana.analyze_eq(solver=select_solver)
        pH = ana.get_pH()
        symbol_of_Ba = ana.__dict__['logbook']['Ba2+']
        self.assertAlmostEqual(ana.__dict__['cation_reference'][symbol_of_Ba]['value'], 1.29 * 10 ** -4, 2)


if __name__ == '__main__':
    unittest.main()
