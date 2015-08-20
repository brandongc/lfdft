import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

class WaveFunctions:
    """Grid and Lagrage function basis wavefunctions"""
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

        for i in xrange(self.fn.nbands):
            s = np.linalg.norm(self.psi[:,i]) * np.sqrt(self.grid.vol)
            self.psi[:,i] /= s

    def iterate(self, hamiltonian):
        A = LinearOperator((self.grid.n, self.grid.n),
                           matvec=hamiltonian.apply,
                           dtype=np.float)
        e_n, self.psi = eigsh(A, k=self.fn.nbands, which='SA',
                              v0=self.psi[:,0])
        self.energy = np.dot(self.fn.f, e_n) * self.grid.vol
        self.orthonormalize()
                                                    
    def calculate_density(self):
        """Calculate a new density from the current wavefunctions"""
        rho = np.zeros(self.grid.n)
        for i in xrange(self.fn.nbands):
            rho[:] += self.fn[i] * self.psi[:,i]**2
        return rho
