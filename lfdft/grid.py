import numpy as np
from scipy.spatial.distance import cdist
from lfdft.io import separator_string

def _print_grid(grid, out):
    out.write(separator_string('Grid'))
    out.write('gpts:    {} {} {}\n'.format(*grid.gpts))
    out.write('spacing: {} {} {}\n'.format(*grid.h))
    out.write('vol:     {}\n'.format(grid.vol))
    out.write('N pts:   {}\n'.format(grid.n))
    

class GridDesc:
    """Descriptor of computational grid.

    todo: document class attributes
    """
    def __init__(self, p, atoms, out):
        self.gpts = (atoms.unitcell / p['grid_spacing']).astype(np.int) + 1
        self.n = self.gpts.prod()
        
        self.h = atoms.unitcell / self.gpts
        self.vol = self.h.prod()

        x = [np.linspace(0, L, N) for L,N in zip(atoms.unitcell, self.gpts)]
        
        X,Y,Z = np.meshgrid(*x, indexing='ij')
        I,J,K = np.meshgrid(np.arange(self.gpts[0]),np.arange(self.gpts[1]),
                            np.arange(self.gpts[2]), indexing='ij')
        self.i = np.vstack((I.ravel(),J.ravel(),K.ravel())).T
        self.p = np.vstack((X.ravel(),Y.ravel(),Z.ravel())).T

        _print_grid(self, out)
        
    def local(self, position, radius):
        r = cdist(self.p, position.reshape(1,-1)).ravel()
        mask = r < radius
        return r[mask], np.arange(self.n)[mask]

