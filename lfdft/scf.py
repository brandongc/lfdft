import numpy as np
from lfdft.io import separator_string

def _print_scf_info(scf, out):
    out.write(separator_string('SCF'))
    out.write('maxiter: {}\n'.format(scf.maxiter))
    out.write('density_error_tol (e / # electrons): {}\n'.format(
        scf.density_error_tol))
    out.write('energy_error_tol (eV / # electrons): {}\n'.format(
        scf.energy_error_tol))

    
class SCF:
    def __init__(self, p, out):
        self.maxiter = p['max_scf_iter']
        self.density_error_tol = p['density_error_tol']
        self.energy_error_tol = p['energy_error_tol']

        self.out = out

        self.total_energy = list()
        self.density_error = 100.0

        _print_scf_info(self, out)
        

    def run(self, wfs, hamiltonian, density, fn):
        """
        1)
        update wfs based on hamiltonian
        orthogonalizes if needed           
        2) 
        calculate new density from updated wfs
        mix with previous density            
        
        3) update the hamiltonian based on the new density
        4) update energies and check convergence
        """
        self.out.write(separator_string('SCF LOOP'))
        self.out.write("{} | {} | {} | {} | {}\n".format(
            'iteration', 'Q', 'Energy', 'Density Error', 'Energy Error'))
        for i in xrange(1, self.maxiter+1):
            wfs.iterate(hamiltonian)
            density.update(wfs)
            hamiltonian.update(density)

            if self.converged(i, hamiltonian, wfs, density):
                break
        self.out.write(separator_string())

    def converged(self, i, hamiltonian, wfs, density):
        energy = (hamiltonian.hartree.e + hamiltonian.xc.e + wfs.energy)
        self.total_energy.append(energy)
        q = density.get_total_charge()
        self.density_error = density.mixer.get_error() / q

        if i >= 3:
            self.energy_error = np.ptp(self.total_energy[-3:]) / q

            txt = "{:3d}  {:0.3f}  {:0.5f}  {:0.5e}  {:0.5e}\n".format(
                i, q, energy, self.density_error, self.energy_error)
            self.out.write(txt)

            if (self.density_error < self.density_error_tol and 
                self.energy_error < self.energy_error_tol):
                converged = True
            else:
                converged = False
        else:
            txt = "{:3d}  {:0.3f}  {:0.5f}  {:0.5e}  {}\n".format(
                i,q,self.total_energy[-1],self.density_error,'   --')
            self.out.write(txt)
            converged = False
            
        return converged
        
        
        
