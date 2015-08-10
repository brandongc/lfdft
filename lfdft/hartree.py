import numpy as np

from lfdft.poisson import Psolver
from lfdft.pseudopotential import calculate_background
from lfdft.units import e2

class Hartree:
    def __init__(self, grid, atoms):
        self.grid = grid
        self.psolver = Psolver(grid)
        self.rho_bg, self.v_bg = calculate_background(grid, atoms)

    def update(self, density):
        b = -4*np.pi * (density.rho - self.rho_bg)
        self.v = e2*self.psolver.solve(b)
        self.e = -0.5*(density.rho * (self.v + self.v_bg)).sum()*self.grid.vol
        
