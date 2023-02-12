from ionlab.IonicEquilibrium.analyze_solution import Lab
from ionlab.IonicEquilibrium.db.matter import Solution, Compound
from ionlab.IonicEquilibrium.db.common_chemicals import *
import ionlab.conf
ionlab.conf.Paths.WOLFRAM_KERNEL = ''
sol = Solution(1)
HI = Compound(cation='H+', anion='I-', cation_count=1, anion_count=1)
# NaOH = Compound(cation='Na+', anion='OH-', cation_count=1, anion_count=1)
sol.add(HI, init_M=0.0065)  # [H3PO4]0 = 0.02(M)
# sol.add(NaOH, init_M=0.1)  # [NaOH]0 = 0.001(M)
lab = Lab(sol)
eq = lab.analyze_eq(solver='gekko')
print(eq)
pH = lab.get_pH()
print(lab.logbook)