import numpy as np
from scipy.sparse.linalg import LinearOperator, lobpcg

class WaveFunctions:
    """Grid and Lagrage function basis wavefunctions
    self.psi.shape = (num_grid_points, num_bands)
    this choice means the current implimentation of
    calculate_density
    and the normalization part of orthonormalize
    have bad memory access patterns.
    This is preferable to doing a copy/transpose of psi for the eigensolver
    """
    def __init__(self, grid, fn):
        self.grid = grid
        self.fn = fn

        self.psi = -0.5 + np.random.random((grid.n, fn.nbands))
        self.orthonormalize()

    def orthonormalize(self):
        """Othogonlize via cholesky factorization and normalize psi"""
        S = np.dot(self.psi.T, self.psi)
        L = np.linalg.cholesky(S)
        Li = np.linalg.inv(L)
        self.psi[:] = np.dot(self.psi, Li)

        # todo: vectorize (bad memory access in current form)
        for i in xrange(self.fn.nbands):
            s = np.linalg.norm(self.psi[:,i]) * np.sqrt(self.grid.vol)
            self.psi[:,i] /= s

    def iterate(self, hamiltonian, tol=2e-8):
        A = LinearOperator((self.grid.n, self.grid.n),
                           matvec=hamiltonian.apply,
                           dtype=np.float)
        e_n, self.psi = lobpcg(A, self.psi, largest=False, tol=tol)
        self.energy = np.dot(self.fn.f, e_n)
        self.orthonormalize()
        #print self.e_n
        #import sys
        #sys.exit()
                        
                            
    def calculate_density(self):
        """Calculate a new density from the current wavefunctions

        Better for the memory access to be inefficient here instead
        of trying to store a transpose of psi after the wavefunctions are 
        iterated.
        """
        rho = np.zeros(self.grid.n)
        for i in xrange(self.fn.nbands):
            rho[:] += self.fn[i] * self.psi[:,i]**2
        return rho
