import numpy as np

from lfdft.pseudopotential import Pseudopotential
from lfdft.hartree import Hartree
from lfdft.operators import Kinetic
from lfdft.xc import XC

class Hamiltonian:
    def __init__(self, grid, setups, atoms, density):

        self.pp = Pseudopotential(grid, atoms, setups)
        self.hartree = Hartree(grid, atoms)
        self.kinetic = Kinetic(grid)
        self.xc = XC(grid)
        self.update(density)

    def update(self, density):
        self.hartree.update(density)
        self.xc.update(density)

    def apply(self, psi_i):
        """matvec routine called by lobpcg eigensolver
        -needs some slicing to keep broadcast rules happy
        -todo: add 'matmat' routine
        """
        return (self.kinetic.apply(psi_i[:,0]) +
                self.pp.apply(psi_i[:,0]) +
                (self.hartree.v + self.xc.v)*psi_i[:,0])
