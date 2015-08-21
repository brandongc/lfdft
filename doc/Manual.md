# lfdft manual

## Options
Options may be set via "key = value" pairs in an input file or as keyword arguments to the LFDFT calculator object.

  * #####max_scf_iter
  number of self consistent field iterations for ground state
   * type: int
   * default: 300

  * #####grid_spacing
  upper limit on distance between grid points in Angstroms
   * type: float
   * default: 0.2

  * #####atomsfile
  name of file with atomic positions (xyz format)
   * type: str
   * default: 'atoms.xyz'

  * #####density_error_tol
  convergence criterial for ground state (difference per electron)
   * type: float
   * default: 1.0e-4

  * #####energy_error_tol
  convergence criterial for ground state (eV difference per electron)
    * type: float
    * default: 5.0e-4

 * #####kinetic_op
 type of kinetic operator. possible options are 'fd' for finite difference and 'lf' for lagrange function
  * type: str
  * default: fd

 * #####txt
 where to write output information. either a file name or 'std' for standard output
  * type: str
  * default: 'std'

 * #####save_psi
 write converged wavefunction to psi.npy
  * type: bool
  * default: False


  
