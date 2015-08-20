import numpy as np

from lfdft.pseudopotential import Pseudopotential
from lfdft.hartree import Hartree
from lfdft.operators import Kinetic, LFKinetic
from lfdft.xc import XC

class Hamiltonian:
    def __init__(self, p, grid, setups, atoms, density):

        # maybe this isn't the best place for this decision???
        if p['kinetic_op'] == 'fd':
            self.kinetic = Kinetic(grid)
        elif p['kinetic_op'] == 'lf':
            self.kinetic = LFKinetic(grid)
        else:
            raise NotImplementedError('incorrect kinetic type: ' +
                                      p['kinetic_op'])
        
        self.pp = Pseudopotential(grid, atoms, setups)
        self.hartree = Hartree(grid, atoms)
        self.xc = XC(grid)
        self.update(density)

    def update(self, density):
        self.hartree.update(density)
        self.xc.update(density)

    def apply(self, psi_i):
        """matvec routine called by eigensolver"""
        return (self.kinetic.apply(psi_i) +
                self.pp.apply(psi_i) +
                (self.hartree.v + self.xc.v)*psi_i)
