import sys

import numpy as np

from lfdft.grid import GridDesc
from lfdft.setups import Setups
from lfdft.occupations import Occupations
from lfdft.density import Density
from lfdft.hamiltonian import Hamiltonian
from lfdft.wfs import WaveFunctions
from lfdft.scf import SCF
from lfdft.atoms import Atoms
from lfdft.io import parse_options, print_parameters, print_atoms

class LFDFT:
    """Main object for a groundstate calculation."""
    def __init__(self, **kwargs):

        p = parse_options(kwargs)
        self.p = p
        if p['txt'] == 'std':
            self.out = sys.stdout
        else:
            self.out = open(p['txt'], 'w', 1)

        print_parameters(p, self.out)
            
        atoms = Atoms(p['atomsfile'])
        print_atoms(atoms, self.out)
        
        self.scf = SCF(p, self.out)
        self.grid = GridDesc(p, atoms, self.out)
        self.setups = Setups(atoms, self.out)
        atoms.set_charges(self.setups)
        self.fn = Occupations(atoms, self.setups, self.out)
        self.density = Density(self.grid, self.setups, atoms)            
        self.hamiltonian = Hamiltonian(p, self.grid, self.setups, atoms,
                                       self.density)
        self.wfs = WaveFunctions(self.grid, self.fn)            
        
    def solve_groundstate_scf(self):
        self.scf.run(self.wfs, self.hamiltonian, self.density, self.fn)

        self.out.write('Hartree energy (eV): {}\n'.format(
            self.hamiltonian.hartree.e))
        self.out.write('XC      energy (eV): {}\n'.format(
            self.hamiltonian.xc.e))
        self.out.write('KS      energy (eV): {}\n'.format(
            self.wfs.energy))
        self.out.write('Total   energy (eV): {}\n'.format(
            self.scf.total_energy[-1]))

        if self.p['save_psi']:
            np.save('psi.npy', self.wfs.psi)
        
        self.out.close()
        return self.scf.total_energy[-1]


