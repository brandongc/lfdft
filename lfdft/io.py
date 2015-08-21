

_default_parameters = {
    'max_scf_iter' : 300,
    'grid_spacing' : 0.2,  #Angstrom
    'atomsfile' : 'atoms.xyz',
    'density_error_tol' : 1.0e-4, # e / electron
    'energy_error_tol' : 5.0e-4,   # eV / electron
    'kinetic_op' : 'fd',  # type of kinetic operator ( 'fd' or 'lf' )
    'txt' : 'std',
    'save_psi' : False
}

def separator_string(txt=None, n=80, symbol='*', pad=1):
    """Print an len(n) string of symbol with optional text in middle"""
    if txt is None:
        return n*symbol + '\n'
    else:
        assert( len(txt) < n-4 )
        n1 = n / 2 - len(txt) / 2
        n2 = n - n1 - len(txt) - 2*pad                
        return n1*symbol + ' ' + txt + ' ' + n2*symbol + '\n'
        
def parse_input_file(infile):
    """Parse key = value pairs from infile.

    -Lines begining with '#' are comments
    -keys may only be set once in input file
    """
    keys_from_file = set()
    options = dict()
    
    for line in infile:
        if '=' in line and line.lstrip()[0] != '#':            
            key, value = line.rstrip().split('=')
            key = key.lower().strip()
            value = value.strip()

            if key in keys_from_file:
                raise ValueError('Duplicate key set in input file: '+key)
            else:
                keys_from_file.add(key)
                options[key] = value
    return options

def parse_options(opts):
    """Get calculation options from 'opts' dict()
    
    -only options in _default_parameters will be recognized    
    -Value's type determined from _default_parameters
    """
    options = dict(_default_parameters)
    for key, value in opts.iteritems():
        if key in _default_parameters:
            default_value = _default_parameters[key]
            options[key] = type(default_value)(value)
    return options                

def print_parameters(parameters, out):
    out.write(separator_string('Parameters'))
    for k, v in parameters.iteritems():
        out.write('{} = {}\n'.format(k,v))

def print_atoms(atoms, out):
    out.write(separator_string('Atoms'))
    out.write('unit cell: {} {} {}\n'.format(*atoms.unitcell))
    for i, (s, p) in enumerate(zip(atoms.symbols, atoms.positions)):
        out.write('{0:4}. {1:2} {2:20.15f} {3:20.15f} {4:20.15f}' \
                  '\n'.format(i+1,s,*p))
