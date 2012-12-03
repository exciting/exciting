import os
import numpy as np
'''
input validation module

provides functions to validate all input variables to the Jacapo calculator.
'''

###########################################3
### Helper functions
##########################################
def get_dacapopath():
    """Return the value of the DACAPOPATH environment variable,
    or, if DACAPOPATH not set, return /usr/share/dacapo-psp/.
    Note that DACAPOPATH must be set for the Fortran code!  """
    import os
    return os.environ.get('DACAPOPATH', '/usr/share/dacapo-psp')

###########################################3
### Validation functions
##########################################
def valid_int(x):
    if (isinstance(x, int) or isinstance(x,np.int32)):
        return True

def valid_float(x):
    return isinstance(x, float)

def valid_int_or_float(x):
    return ((isinstance(x, int) or isinstance(x,np.int32))
            or isinstance(x, float))

def valid_boolean(x):
    return isinstance(x, bool)

def valid_str(x):
    return isinstance(x, str)

def valid_atoms(x):
    import ase
    return isinstance(x, ase.Atoms)

def valid_pw(x):
    return (valid_int_or_float(x) and x>0 and x<2000)

def valid_dw(x):
    return (valid_int_or_float(x) and x>0 and x<2000)

def valid_xc(x):
    return (x in ['PW91', 'PBE', 'revPBE', 'RPBE', 'VWN'])

def valid_nbands(x):
    return valid_int(x)

def valid_ft(x):
    return(valid_float, x)

def valid_spinpol(x):
    return valid_boolean(x)

def valid_fixmagmom(x):
    return valid_float(x)

def valid_symmetry(x):
    return valid_boolean(x)

def valid_calculate_stress(x):
    return valid_boolean(x)

def valid_kpts(x):
    if isinstance(x, str):
        return x in ['cc6_1x1',
                     'cc12_2x3',
                     'cc18_sq3xsq3',
                     'cc18_1x1',
                     'cc54_sq3xsq3',
                     'cc54_1x1',
                     'cc162_sq3xsq3',
                     'cc162_1x1']
    x = np.array(x)
    #empty arg is no good
    if x.shape == ():
        return False
    #monkhorst-pack
    elif x.shape == (3,) and ((x.dtype == 'int32') or (x.dtype == 'int64')):
        return True
    #user-defined list
    elif x.shape[1] == 3 and (str(x.dtype))[0:7] == 'float64':
        return True
    else:
        return False

def valid_dipole(x):
    if valid_boolean(x):
        return True
    #dictionary passed in. we need to check the keys
    valid_keys = {'status':valid_boolean,
                  'mixpar':valid_float,
                  'initval':valid_float,
                  'adddipfield':valid_float,
                  'position':valid_float}
    for key in x:
        if key not in valid_keys:
            return False
        else:
            if x[key] is not None:
                if not valid_keys[key](x[key]):
                    return False
    return True

def valid_nc(x):
    #todo check for read/write access?
    return valid_str(x)

def valid_status(x):
    return valid_str(x)

def valid_pseudopotentials(x):
    #todo check that keys are symbols or numbers
    #todo check that psp files exist

    return True

    dacapopath = get_dacapopath()
    if dacapopath is None:
        raise Exception, 'No $DACAPOPATH found. please set it in .cshrc or .bashrc'

    from ase.data import chemical_symbols
    for key in x:
        if valid_str(key):
            if key not in chemical_symbols:
                return False
        elif not (valid_int(key) and key > 0 and key < 112):
            return False

        #now check for existence of psp files
        psp = x[key]
        if not (os.path.exists(psp)
                or os.path.exists(os.path.join(dacapopath, psp))):
            return False
    return True

def valid_extracharge(x):
    return valid_float(x)

def valid_extpot(x):
    grids = get_fftgrid()
    if (x.shape == np.array(grids['soft'])).all():
        return True
    else:
        return False

def valid_ascii_debug(x):
    return (x.strip() in ['Off', 'MediumLevel', 'HighLevel'])

def valid_ncoutput(x):
    if x is None:
        return
    valid_keys = ['wf', 'cd', 'efp', 'esp']

    for key in x:
        if key not in valid_keys:
            return False
        else:
            if x[key] not in ['Yes', 'No']:
                return False
    return True

def valid_ados(x):
    if x is None:
        return
    valid_keys = ['energywindow',
                  'energywidth',
                  'npoints',
                  'cutoff']
    for key in x:
        if key not in valid_keys:
            print '%s not in %s' % (key, str(valid_keys))
            return False
        if key == 'energywindow':
            if not len(x['energywindow']) == 2:
                print '%s is bad' % key
                return False
        if key == 'energywidth':
            if not valid_float(x['energywidth']):
                print key, ' is bad'
                return False
        elif key == 'npoints':
            if not valid_int(x['npoints']):
                print key, ' is bad'
                return False
        elif key == 'cutoff':
            if not valid_float(x['cutoff']):
                print key, ' is bad'
                return False
    return True


def valid_decoupling(x):
    if x is None:
        return
    valid_keys = ['ngaussians', 'ecutoff', 'gausswidth']
    for key in x:
        if key not in valid_keys:
            return False
        elif key == 'ngaussians':
            if not valid_int(x[key]):
                print key
                return False
        elif key == 'ecutoff':
            if not valid_int_or_float(x[key]):
                return False
        elif key == 'gausswidth':
            if not valid_float(x[key]):
                print key, x[key]
                return False
    return True

def valid_external_dipole(x):
    if x is None:
        return
    if valid_float(x):
        return True

    valid_keys = ['value', 'position']

    for key in x:
        if key not in valid_keys:
            return False
        if key == 'value':
            if not valid_float(x['value']):
                return False
        elif key == 'position':
            if not valid_float(x['position']):
                return False
    return True

def valid_stay_alive(x):
    return valid_boolean(x)

def valid_fftgrid(x):
    valid_keys = ['soft', 'hard']
    for key in x:
        if key not in valid_keys:
            return False
        if x[key] is None:
            continue

        grid = np.array(x[key])
        if (grid.shape != (3,) and grid.dtype != 'int32'):
            return False
    return True

def valid_convergence(x):
    valid_keys = ['energy',
                  'density',
                  'occupation',
                  'maxsteps',
                  'maxtime']
    for key in x:
        if key not in valid_keys:
            return False
        if x[key] is None:
            continue
        if key == 'energy':
            if not valid_float(x[key]):
                return False
        elif key == 'density':
            if not valid_float(x[key]):
                return False
        elif key == 'occupation':
            if not valid_float(x[key]):
                return False
        elif key == 'maxsteps':
            if not valid_int(x[key]):
                return False
        elif key == 'maxtime':
            if not valid_int(x[key]):
                return False
    return True

def valid_charge_mixing(x):
    valid_keys = ['method',
                  'mixinghistory',
                  'mixingcoeff',
                  'precondition',
                  'updatecharge']

    for key in x:
        if key not in valid_keys:
            return False
        elif key == 'method':
            if x[key] not in ['Pulay']:
                return False
        elif key == 'mixinghistory':
            if not valid_int(x[key]):
                return False
        elif key == 'mixingcoeff':
            if not valid_float(x[key]):
                return False
        elif key == 'precondition':
            if x[key] not in ['Yes', 'No']:
                return False
        elif key == 'updatecharge':
            if x[key] not in ['Yes', 'No']:
                return False
    return True

def valid_electronic_minimization(x):
    valid_keys = ['method', 'diagsperband']
    for key in x:
        if key not in valid_keys:
            return False
        elif key == 'method':
            if x[key] not in ['resmin',
                              'eigsolve',
                              'rmm-diis']:
                return False
        elif key == 'diagsperband':
            if not (valid_int(x[key]) or x[key] is None):
                return False
    return True

def valid_occupationstatistics(x):
    return (x in ['FermiDirac', 'MethfesselPaxton'])


def valid_mdos(x):
    return True

def valid_psp(x):

    dacapopath = get_dacapopath()
    
    valid_keys = ['sym','psp']
    if x is None:
        return True
    for key in x:
        if key not in valid_keys:
            return False
        if not valid_str(x[key]):
            return False
        if key == 'sym':
            from ase.data import chemical_symbols
            if key not in chemical_symbols:
                return False
        if key == 'psp':
            
            if os.path.exists(x['psp']):
                return True

            if os.path.exists(os.path.join(dacapopath, x['psp'])):
                return True
            #psp not found
            return False
