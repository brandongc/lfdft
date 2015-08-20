from lfdft import LFDFT

calc = LFDFT(atomsfile='li15.xyz',
             grid_spacing=0.3,
             txt='li15.txt')
total_energy = calc.solve_groundstate_scf()
