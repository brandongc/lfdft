from lfdft import LFDFT

calc = LFDFT(atomsfile='c60.xyz',
             grid_spacing=0.3,
             txt='c60.txt')
total_energy = calc.solve_groundstate_scf()
