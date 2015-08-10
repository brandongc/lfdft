from __future__ import print_function
import argparse
import sys

from lfdft.atoms import Atoms
from lfdft.io import parse_input, print_parameters, print_atoms
from lfdft.calculator import Calculator

_description = "Density functional theory with Lagrange function basis"

def main():
    """Parse command line options for DFT calculation."""
    parser = argparse.ArgumentParser(description=_description)
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout)         
    args = parser.parse_args()
    
    p = parse_input(args.infile)
    args.infile.close()
    out = args.outfile
    
    print_parameters(p, out)

    atoms = Atoms(p['structure_file'])
    print_atoms(atoms, out)

    calc = Calculator(p, atoms, out)
    calc.solve_scf()
    
    out.close()
