# setup the geometry
# must manually add cell information to the resulting xyz file (2nd line)
import ase
from ase.cluster.cubic import BodyCenteredCubic

surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
layers = [3, 3, 3]
lc = 3.49000
atoms = BodyCenteredCubic('Li', surfaces, layers, latticeconstant=lc)

atoms.set_cell((20,20,20))
atoms.center()
print atoms.cell
ase.io.write('li15.xyz', atoms)
