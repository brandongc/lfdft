import argparse
from lfdft.io import parse_input_file
from lfdft.calculator import LFDFT

_description = "Density functional theory with Lagrange function basis"

def main():
    """Parse command line options for DFT calculation."""
    parser = argparse.ArgumentParser(description=_description)
    parser.add_argument('infile', type=argparse.FileType('r'),
                        help='input file for calculation (required)')
    args = parser.parse_args()

    options = parse_input_file(args.infile)    
    args.infile.close()
    
    calc = LFDFT(**options)
    calc.solve_groundstate_scf()
