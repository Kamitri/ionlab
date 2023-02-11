# Get Started

Ionlab is a **very** simple python package to determine the composition at equilibrium for a solution containing different species. Ionlab assumes the **activity coeffecient of all ions to be 1**, and all other factors considered ideal. Complexation and precipitation are also taken into account, however computation may be inaccurate or may not be possible at all. 
Equation solving is handled by **GEKKO (pure python)**, or Wolfram Mathematica.
To get started, download ionlab through the Python package manager:

    pip install ionlab
Carry out a simple experiment:

	import ionlab as il
    H3PO4 = il.Compound(cation='H+', anion='PO43-', cation_count=3, anion_count=1)
    sol = il.Solution(1)  # Create a solution with a volume of 1L
    sol.add(H3PO4, init_M=0.2)  # [H3PO4]0 = 0.2(M)
    lab = il.Lab(sol)  
    eq = lab.analyze_eq(solver='gekko')  
    print(eq)
    >>> {'xani0st0': 0.16610105407, 'h': 0.033899007588, 'pH': 1.46981301581946}
    
A more complicated example for the system consisting of satured BaCO<sub>3</sub>

    sol = Solution(1)  
    sol.add(BaCl2, init_M=1)  
    sol.add(Na2CO3, init_M=1)  
    lab = Lab(sol)  
    print(lab.analyze_eq(solver='wolfram'))  # GEKKO would fail here
    >>> {'h': 1.099744438343427e-10, 'xani1st0': 2.238882341544366e-08, 'xcat0st0': 0.0001295842673586922, 'xcat1st0': 2.0, 'xcat0ani1sl': 0.999870364297433, 'xani0st0': 2.199488876686854e-17, 'pH': 9.958708225671051}
    pH = lab.get_pH()
    print(pH)
    >>> 9.958708225671051
To see what the symbols mean, please use the following:

    print(lab.logbook)
    >>> {'H+': 'h', 'Cl-': 'xani0st1', 'HCl': 'xani0st0', 'CO32-': 'xani1st2', 'H2CO3': 'xani1st0', 'HCO3-': 'xani1st1', 'Ba2+': 'xcat0st0', 'BaOH+': 'xcat0st1', 'Na+': 'xcat1st0', 'xcat0ani1sl': 'BaCO3'}

## Switching solvers
I strongly recommend using **Wolfram Mathematica** instead of GEKKO as it is much faster and more reliable at solving complex systems (also because I am a horrible mathematician and didn't properly handle equation solving). If you have Mathematica installed, do the following step to configure the wolfram kernel and utilize it instead.

    import conf
    conf.Paths.WOLFRAM_KERNEL = '<.../WolframKernel.exe>'  # Path to Wolfram kernel goes here
    ...
    lab.analyze_eq(solver='wolfram')  # Use wolfram to solve equation
    
## More compounds
Some common compounds have been provided for quick access. Type `cc.` to see a list of common chemicals. Compounds can be declared like below (*untuitive and inconsistent, I know*):

    NH3 = Compound(formula='NH3')
    CdBr2 = Compound(cation='Cd2+', anion='Br-', cation_count=1, anion_count=2)
    H2SO4 = Compound(formula='H2SO4', cation_count=2, anion_count=1)
    H2SO4 = Compound(cation='H+', anion='SO42-', cation_count=2, anion_count=1)
    
# Limitations
I made this package when I wasn't the greatest programmer, nor the best chemist. This still hold true, but at least now I can ask someone, maybe. Please be aware that this package may be **buggy**, **incomplete**, **inconsistent**. Honestly, thanks for the reading the docs but trying out this package may **not be the best experience**. But well, if you want to explore, go for it. It is unknown whether I will ever update this package. If maybe someone need it, perhaps.
**Current limitations**
 - Computation may fail for systems that are too complex, or... just fail in general.
 - Computation may be mildy or wildy inaccurate for systems involving complexation/precipitation.
 - Unintuitive design for adding/adjusting chemicals (data).
 - Changing chemicals' constant values during runtime is not available.
 - The ideal assumption is utilized.
 - GEKKO solver (pure python) may be slower at calculating and easier to fail for complex systems.