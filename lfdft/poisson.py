import numpy as np
from scipy.sparse.linalg import bicgstab, LinearOperator
from lfdft.operators import Laplacian

class Psolver:
    def __init__(self, grid, tol=1e-8):
        self.tol = tol
        self.x = np.zeros(grid.n)
        nabla = Laplacian(grid)
        self.A = LinearOperator((grid.n, grid.n),
                                matvec = nabla.apply,
                                dtype=np.float)
        
    def solve(self, b):
        self.x, info = bicgstab(self.A,b,x0=self.x,tol=self.tol)        
        return self.x
