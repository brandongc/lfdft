import pkg_resources
import itertools
import numpy as np
from string import maketrans
from scipy.interpolate import UnivariateSpline as spline
from scipy.integrate import simps
from units import e2

from lfdft.io import separator_string
from lfdft.ylm import ylm

def _print_pseudopotential(setup, out):
    out.write('symbol: {}\n'.format(setup.symbol))
    out.write('dx:     {}\n'.format(setup.dx))
    out.write('charge: {}\n'.format(setup.q))
    out.write('radius: {}\n'.format(setup.radius))
    out.write('lmax:   {}\n\n'.format(setup.lmax))
    # for debuging
    #out.write('uVu:\n')
    #for lm, uVu in zip(setup.LM, setup.uVu):
    #    out.write('{} {}\n'.format(lm,uVu))
    
class Setup:
    def __init__(self, symbol):
        """Atomic Setup (Pseudopotential data)
        data read from data/*_pp.dat files.
        TODO: expand documentation & cleanup
        """
        self.lref = 1

        fname = 'data/' + symbol.lower() + '_pp.dat'                
        datafile = pkg_resources.resource_stream(__name__, fname)
        self.symbol = symbol
        self.fname = fname
        
        # python doesn't like 1D-2 formating (convert to 1E-2)
        t = maketrans('D', 'E')
        data = [l.translate(t).split() for l in datafile]
        
        line=0
        self.n = int(data[line][0])
        self.dx = float(data[line][1])
        self.lmax = int(data[line][2])
        self.q = float(data[line][3])
        
        line=1 
        self.rl = np.array([float(r) for r in data[line]]) 
        self.radius = max(self.rl)
        
        m = self.n + 3
        v = np.array([[float(x) for x in line] for line in data[2:m]]).T
        u = np.array([[float(x) for x in line] for line in data[m::]]).T
        
        x = np.arange(self.n+1) * self.dx

        # create splines
        self.v_spline = {l:spline(x, v[l+1], s=0)
                         for l in xrange(self.lmax+1)}
        self.u_spline = {l:spline(x, u[l+1], s=0)
                         for l in xrange(self.lmax+1)}

        self.vnl_spline = {}
        self.psi_spline = {}

        psi = np.zeros_like(x)
        vnl = np.zeros_like(x)
            
        for l in xrange(self.lmax+1):
            psi[1::] = u[l+1,1::]/(x[1::]**(l+1)*np.sqrt((2*l+1)/(4*np.pi)))
            vnl[1::] = (v[l+1,1::]-v[self.lref+1, 1::])*psi[1::]
            psi[0] = psi[1]
            vnl[0] = vnl[1]

            self.psi_spline[l] = spline(x,psi,s=0)
            self.vnl_spline[l] = spline(x,vnl,s=0)
            

        self.L = [l for l in xrange(self.lmax+1) if l != self.lref]        
        self.LM = list(itertools.chain.from_iterable(
                range(l**2, (l+1)**2) for l in self.L))
        self.uVu = self.calculate_uVu()
        
    def calculate_uVu(self, npoints=1000):
        """Redundancy is used here to allow for numpy array operations instead
        of python loops in application of non-local pseudopotential.

        A typical radius is < 1 Angstrom so 1000 points should always be 
        sufficiently accurate
        """
        x = np.linspace(0, self.radius, npoints)
        uVu = np.zeros(len(self.LM))

        lm_index = 0
        for l in self.L:
            uVuL = 4*np.pi*simps(self.psi(l,x)*self.vnl(l,x)*x**2 / (2*l+1), x)
            for lm in xrange(l**2, (l+1)**2):
                uVu[lm_index] = uVuL
                lm_index += 1
        return uVu

    def calculate_uV(self, x, r):
        uV = np.zeros((len(self.LM),len(x)))
        for l in self.L:
            uVrl = self.vnl(l, r)
            for i, lm in enumerate(self.LM):
                uV[i,:] = uVrl * ylm(x,lm+1)
        return uV
                  
    def u(self, l, x):
        m = x < self.radius
        u = np.zeros_like(x)
        u[m] = self.u_spline[l]( x[m] )
        return u
    
    def v(self, x, l=None):
        "v(r) = -q/r for r > R"
        if l is None:
            l = self.lref
        m = x < self.radius
        v = np.empty_like(x)
        v[m] = self.v_spline[l](x[m])
        m = np.invert(m)
        v[m] = -self.q * e2 / x[m]
        return v
        
    def vnl(self, l, x):
        return self.vnl_spline[l](x)

    def psi(self, l, x):
        return self.psi_spline[l](x)



def Setups(atoms, out):
    """Returns dict of {symbol : Setup(symbol)}."""
    out.write(separator_string('Pseudopotentials'))

    setups = dict()
    for s in set(atoms.symbols):
        setups[s] = Setup(s)
        _print_pseudopotential(setups[s], out)

    return setups
    

