import numpy as np
from scipy.spatial.distance import cdist

from lfdft.mixer import Mixer

def _gaussian_density(grid, atoms, beta=1.5):
    """Generate a density on grid using gaussians centered on the atoms.
    Gaussians are scaled by the charge on the atom to improve the quality
    of the initial guess.
    """    
    rho = np.zeros(grid.n)
    for pos, q in zip(atoms.positions, atoms.charges):
        r = cdist(grid.p, pos.reshape(1,-1)).ravel()
        rho += q * np.exp(-beta*r**2)
    rho /= rho.sum() * grid.vol / atoms.total_electrons
    return rho

class Density:
    def __init__(self, grid, setups, atoms):
        self.grid = grid
        self.rho = _gaussian_density(grid, atoms)
        self.mixer = Mixer(grid)
        self.mixer.mix(self.rho)
        
    def update(self, wfs):
        self.rho = wfs.calculate_density()
        self.mixer.mix(self.rho)

    def get_total_charge(self):
        return self.rho.sum() * self.grid.vol
