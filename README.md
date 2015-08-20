# lfdft
###Density Functional Theory with Finite Difference and Lagrange Function Basis Set
A easy to run, understand and modify package allowing rapid implementation and testing of computational methods.

## Current features
  * LDA approximation for exchange and correlation
  * non-periodic systems only
  * finite difference or Lagrange function kinetic operator
  * built-in set of Pseudopotentials:
    * Ag, Al, Au, Be, B, Br, Cl, C, Cs, Cu, F, Ge, He, H, K, Li, Mg, Na, Ni, N, O, Pd, P, Rb, Si, S, Ti, V, W

## Coming soon
In general these features have been implemented, but need to be integrated with this code base and/or tested.
  * Real-time time-dependent density functional theory
   * Time-dependent current (including non-local contribution)
   * Erhenfest dynamics
   * Laser Pulse
  * Periodic boundary conditions
  * Atomic Orbitals (LCAO)
  * Atomic forces (e.g. geometry optimization)
  * Complex potential transport

## Install
Requirements: NumPy and SciPy.
Installing in a [virtualenv](https://virtualenv.pypa.io) is recommended.

## Usage
After e.g. `python setup.py develop` to install the package run from the command line with `$ lfdft input_file`

Atomic geometry is specified by [XYZ format](https://en.wikipedia.org/wiki/XYZ_file_format) where the comment line contains the unit cell size in Angstroms.
Example C2H4 in 10.0 x 10.0 x 10.0 Angstrom box:

    6
    10.0 10.0 10.0
    C       5.000000000000000      5.000000000000000      5.667480000000000
    C       5.000000000000000      5.000000000000000      4.332520000000000
    H       5.000000000000000      5.922832000000000      6.237695000000000
    H       5.000000000000000      4.077168000000000      6.237695000000000
    H       5.000000000000000      5.922832000000000      3.762305000000000
    H       5.000000000000000      4.077168000000000      3.762305000000000

Note: the computational grid is defined on [0, L] where L is length of the cell in a given direction. Therefore it is important to ensure the atoms are centered in the compuational box.
