import numpy as np
from scipy.spatial.distance import cdist
from scipy.special import erf
from lfdft.units import e2


def _calculate_local_pseudopotential(grid, atoms, setups):
    v = np.zeros(grid.n)    
    min_r = np.finfo(np.double).tiny
    for p, s in zip(atoms.positions, atoms.symbols):
        r = cdist(grid.p, p.reshape(1,-1)).ravel()
        r[r<min_r] = min_r
        v += setups[s].v(r)
    return v

def calculate_background(grid, atoms, beta=7.0, rmin=1e-16):
    n, v = np.zeros(grid.n), np.zeros(grid.n)

    for p, q in zip(atoms.positions, atoms.charges):
        r = cdist(grid.p, p.reshape(1,-1)).ravel()        
        n += q * np.exp(-beta*r**2) * (beta/np.pi)**1.5
        r[r<rmin] = rmin
        v += erf(np.sqrt(beta)*r) / r*e2*q
    return n, v


class NonLocalPseudopotential:
    def __init__(self, grid, setup, position):
        self.grid = grid
        self.setup = setup
        self.position = position

        r, self.points = grid.local(position, setup.radius)
        x = grid.p[self.points, :] - position
        self.uV = self.setup.calculate_uV(x, r)
        self.uVu = setup.uVu
        
    def apply(self, psi_i, vnl):
        uV_psi = (self.uV * psi_i[self.points]).sum(axis=-1)
        c = self.grid.vol / self.uVu
        vnl[self.points] += (c * uV_psi * self.uV.T).sum(axis=-1)

class Pseudopotential:
    def __init__(self, grid, atoms, setups):

        self.v_local = _calculate_local_pseudopotential(grid, atoms, setups)
        self.v_nonlocal = [NonLocalPseudopotential(grid,setups[s],p) for
                           p,s in zip(atoms.positions, atoms.symbols)]

    def apply(self, psi_i):
        vnl = np.zeros_like(psi_i)
        for nlpp in self.v_nonlocal:
            nlpp.apply(psi_i, vnl)
        return vnl + self.v_local * psi_i
