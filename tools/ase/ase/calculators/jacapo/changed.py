import numpy as np

import logging
log = logging.getLogger('Jacapo')

'''
provides functions to determine if an input parameter has changed.
'''

#######################################################################
#### changed functions

def kpts_changed(calc, x):
    '''
    check if kpt grid has changed.

    we have to take care to generate the right k-points from x if
    needed. if a user provides (4,4,4) we need to generate the MP
    grid, etc...

    Since i changed the MP code in set_kpts, there is some
    incompatibility with old jacapo calculations and their MP
    grids.
    '''
    #chadi-cohen
    if isinstance(x, str):
        exec('from ase.dft.kpoints import %s' % x)
        listofkpts = eval(x)
    #monkhorst-pack grid
    elif np.array(x).shape == (3,):
        from ase.dft.kpoints import monkhorst_pack
        N1, N2, N3 = x
        listofkpts = monkhorst_pack((N1, N2, N3))
    #user-defined list is provided
    elif len(np.array(x).shape) == 2:
        listofkpts = np.array(x)
    else:
        raise Exception, 'apparent invalid setting for kpts'

    grid = calc.get_kpts()
    
    if grid.shape != listofkpts.shape:
        return True

    if (abs(listofkpts - grid) < 1e-6).all():
        return False
    else:
        return True

def electronic_minimization_changed(calc, x):
    myx = calc.get_electronic_minimization()

    for key in myx:
        if myx[key] != x[key]:
            print key, myx[key], ' changed to ', x[key]
            return True
    return False

def spinpol_changed(calc, x):
    if x != calc.get_spinpol():
        return True
    else:
        return False

def symmetry_changed(calc, x):
    if x != calc.get_symmetry():
        return True
    else:
        return False

def xc_changed(calc, x):
    if x != calc.get_xc():
        return True
    return False

def calculate_stress_changed(calc, x):
    if x != calc.get_calculate_stress():
        return True
    return False

def ados_changed(calc, x):
    ados = calc.get_ados()

    #ados may not be defined, and then None is returned
    if ados is None and x is None:
        return False
    elif ados is None and x is not None:
        return True
    elif ados is not None and x is None:
        return True

    #getting here means ados and x are not none so we compare them
    for key in x:
        try:
            if x[key] != ados[key]:
                return True
        except ValueError:
            if (x[key] != ados[key]).all():
                return True
    return False

def convergence_changed(calc, x):
    conv = calc.get_convergence()
    for key in x:
        if x[key] != conv[key]:
            return True
    return False

def charge_mixing_changed(calc, x):
    cm = calc.get_charge_mixing()
    if x is None and cm is None:
        return False
    else:
        return True
        
    for key in x:
        if x[key] != cm[key]:
            return True
    return False

def decoupling_changed(calc, x):
    pars = calc.get_decoupling()
    for key in x:
        if x[key] != pars[key]:
            return True
    return False

def dipole_changed(calc, x):

    pars = calc.get_dipole() #pars stored in calculator

    # pars = False if no dipole variables exist
    if (pars is False and x is False):
        return False #no change
    elif (pars is False and x is not False):
        return True

    # both x and pars is a dictionary
    if (type(pars) == type(dict) and
        type(pars) == type(x)):
        for key in x:
            if key == 'position':    # dipole layer position is never writen to the nc file
                print 'need to do something special'
                continue
            if x[key] != pars[key]:
                return True

    #nothing seems to have changed.
    return False

def extpot_changed(calc, x):
    extpot = calc.get_extpot()
    if (x == extpot).all():
        return False
    return True

def fftgrid_changed(calc, x):
    validkeys = ['soft', 'hard']

    myx = calc.get_fftgrid()
    if (myx['soft'] == x['soft'] and myx['hard'] == x['hard']):
        return False
    else:
        return True


def nbands_changed(calc, x):
    if calc.get_nbands() == x:
        return False
    else:
        return True

def occupationstatistics_changed(calc, x):
    if calc.get_occupationstatistics() == x:
        return False
    else:
        return True

def pw_changed(calc, x):
    if calc.get_pw() == x:
        return False
    else:
        return True

def dw_changed(calc, x):
    if calc.get_dw() == x:
        return False
    else:
        return True

def ft_changed(calc, x):
    if calc.get_ft() == x:
        return False
    else:
        return True
    
def mdos_changed(calc,x):

    myx = calc.get_mdos()

    log.debug('myx = %s' % str(myx))
    log.debug('x = %s' % str(x))

    if x is None and myx is None:
        return False
    elif ((x is None and myx is not None)
        or (x is not None and myx is None)):
        return True
    else:
        for key in x:
            if x[key] != myx[key]:
                return True
    return False

def pseudopotentials_changed(calc,x):

    mypsp = calc.get_pseudopotentials()

    if len(mypsp) != len(x):
        return True

    for key in x:
        if key not in mypsp:
            return True
        if mypsp[key] != x[key]:
            return True

    for key in mypsp:
        if key not in x:
            return True
        if mypsp[key] != x[key]:
            return True
    return False

def status_changed(calc,x):
    if calc.get_status() != x:
        return True
    return False
