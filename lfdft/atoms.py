import numpy as np

class Atoms:
    def __init__(self, atomfile):
        """atomfile should be in .xyz style format.
        1. # of atoms
        2. unit cell (3 floats)
        3. element x y z
        4. ...


        NOTE: line 2 is usually a comment line, but an error will be thrown if
        a unit cell is not provided.
        """        
        with open(atomfile, 'r') as f:
            natoms = int(f.readline().strip())
            unitcell = [float(x) for x in f.readline().split()]
            
            assert len(unitcell) == 3, ('line 2 of {} must contain unit cell' \
                                        '(3 values)'.format(atomfile) )

            positions = []
            symbols = []

            for i in xrange(natoms):
                symbol, x, y, z = f.readline().split()[:4]
                symbols.append(symbol.lower())
                positions.append([float(x), float(y), float(z)])


            self.natoms = natoms
            self.positions = np.array(positions)
            self.symbols = symbols
            self.unitcell = np.array(unitcell)
                
    def set_charges(self, setups):
        self.charges = [setups[s].q for s in self.symbols]
        self.total_electrons = sum(self.charges)

            
