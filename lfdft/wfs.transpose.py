import numpy as np
from scipy.sparse.linalg import LinearOperator, lobpcg

class WaveFunctions:
    def __init__(self, grid, fn):
        self.grid = grid
        self.fn = fn

        self.psi = -0.5 + np.random.random((fn.nbands, grid.n))
        print 'psi', self.psi.shape
        self.orthonormalize()

    def orthonormalize(self):
        """Othogonlize via cholesky factorization and normalize psi"""
        S = np.dot(self.psi, self.psi.T)
        L = np.linalg.cholesky(S)
        Li = np.linalg.inv(L)
        self.psi[:] = np.dot(Li, self.psi)

        for i in xrange(self.fn.nbands):
            s = np.linalg.norm(self.psi[i,:]) * np.sqrt(self.grid.vol)
            self.psi[i,:] /= s

    def iterate(self, hamiltonian, tol=2e-8):
        A = LinearOperator((self.grid.n, self.grid.n),
                           matvec=hamiltonian.apply,
                           dtype=np.float)
        e_n, v = lobpcg(A, self.psi.T, largest=False, tol=tol)
        print v.shape
        self.psi[:] = v.T
        self.e_n = np.dot(self.fn.f, e_n)
        self.orthonormalize()
        print self.e_n
        import sys
        sys.exit()
                        
                            
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
