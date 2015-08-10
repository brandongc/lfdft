import numpy as np
from scipy.ndimage.filters import convolve, convolve1d
from lfdft.units import h2m

# Expansion coefficients for finite difference Laplacian.  The numbers are      
# from J. R. Chelikowsky et al., Phys. Rev. B 50, 11355 (1994):                 
_laplace_coef = [
    [0.],
    [-2.0, 1.0],
    [-5.0/2, 4.0/3, -1.0/12],
    [-49.0/18, 3.0/2, -3.0/20, 1.0/90],
    [-205.0/72, 8.0/5, -1.0/5, 8.0/315, -1.0/560],
    [-5269.0/1800, 5.0/3, -5.0/21, 5.0/126, -5.0/1008, 1.0/3150],
    [-5369.0/1800, 12.0/7, -15.0/56, 10.0/189, -1.0/112, 2.0/1925, -1.0/16632]
]

class Laplacian:
    def __init__(self, grid, order=4):
        self.grid = grid
        l = np.array(_laplace_coef[order])
        s = np.empty(2*order+1)
        s[0:order+1] = l[::-1]
        s[order::] = l
        self.stencil = np.array([s / h**2 for h in grid.h])

    def apply(self, x):
        y = x.reshape(self.grid.gpts)
        A = (convolve1d(y, self.stencil[0], axis=0, mode='constant') +
             convolve1d(y, self.stencil[1], axis=1, mode='constant') +
             convolve1d(y, self.stencil[2], axis=2, mode='constant'))
        return A.reshape(self.grid.n)

class Kinetic(Laplacian):
    """Finite difference kinetic energy operator """
    def __init__(self, grid):
        Laplacian.__init__(self,grid)
        self.stencil *= -h2m

class Metric:
    """Kerker like weighting for density similarity. 
    Weight determines how much to bias away from low-frequency fluctuations.
    todo: add reference
    """
    def __init__(self, grid, weight):
        self.grid = grid

        a = 0.125 * (weight + 7)
        b = 0.0625 * (weight - 1)
        c = 0.03125 * (weight - 1)
        d = 0.015625 * (weight - 1)
                        
        coef = [a,
                b, b, b, b, b, b,
                c, c, c, c, c, c, c, c, c, c, c, c,
                d, d, d, d, d, d, d, d]
        offs = [(0, 0, 0),
                (-1, 0, 0), (1, 0, 0),                 #b
                (0, -1, 0), (0, 1, 0),                 #b
                (0, 0, -1), (0, 0, 1),                 #b
                (1, 1, 0), (1, 0, 1), (0, 1, 1),       #c
                (1, -1, 0), (1, 0, -1), (0, 1, -1),    #c
                (-1, 1, 0), (-1, 0, 1), (0, -1, 1),    #c
                (-1, -1, 0), (-1, 0, -1), (0, -1, -1), #c
                (1, 1, 1), (1, 1, -1), (1, -1, 1),     #d
                (-1, 1, 1), (1, -1, -1), (-1, -1, 1),  #d
                (-1, 1, -1), (-1, -1, -1)              #d
                ]

        s = np.zeros((3,3,3))
        m = 1
        for c, o in zip(coef, offs):
            s[np.array(o)+m] = c
        self.stencil = s


    def apply(self, x):
        return convolve(x.reshape(self.grid.gpts),
                        self.stencil, mode='constant').reshape(self.grid.n)

