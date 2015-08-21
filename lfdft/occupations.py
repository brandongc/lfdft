import numpy as np
from lfdft.io import separator_string

class Occupations:
    def __init__(self, atoms, setups, out):
        qtotal = sum( setups[s].q for s in atoms.symbols )
        if qtotal % 2:
            self.nbands = int(qtotal/2 + 1)
            self.f = 2.0 * np.ones(self.nbands)
            self.f[-1] = 1.0
        else:
            self.nbands = int(qtotal/2)
            self.f = 2.0 * np.ones(self.nbands)


        #out.write(separator_string(''))
        out.write('total electrons:  {}\n'.format(qtotal))
        out.write('electronic bands: {}\n'.format(self.nbands))
                  
                  
    def __getitem__(self, i):
        return self.f[i]
