import numpy as np
from ase import Atoms
from ase.io import write, read
from ase.test import NotAvailable

a = 5.0
d = 1.9
c = a / 2
atoms = Atoms('AuH',
              positions=[(c, c, 0), (c, c, d)],
              cell=(a, a, 2 * d),
              pbc=(0, 0, 1))
extra = np.array([ 2.3, 4.2 ])
atoms.set_array('extra', extra)
atoms *= (1, 1, 2)
images = [atoms.copy(), atoms.copy()]
r = ['xyz', 'traj', 'cube', 'pdb', 'cfg', 'struct', 'cif', 'gen']
try:
    import Scientific
    version = Scientific.__version__.split('.')
    print 'Found ScientificPython version: ',Scientific.__version__
    if map(int,version) < [2,8]:
        print 'ScientificPython 2.8 or greater required for numpy support in NetCDF'
        #raise NotAvailable('ScientificPython version 2.8 or greater is required')
except (ImportError, NotAvailable):
    print 'No Scientific python found. Check your PYTHONPATH'
    #raise NotAvailable('ScientificPython version 2.8 or greater is required')
else:
    r += ['etsf']
w = r + ['xsf', 'findsym']
try:
    import matplotlib
except ImportError:
    pass
else:
    w += ['png', 'eps']

for format in w:
    print format, 'O',
    fname1 = 'io-test.1.' + format
    fname2 = 'io-test.2.' + format
    write(fname1, atoms, format=format)
    if format not in ['cube', 'png', 'eps', 'cfg', 'struct', 'etsf', 'gen']:
        write(fname2, images, format=format)

    if format in r:
        print 'I'
        a1 = read(fname1)
        assert np.all(np.abs(a1.get_positions() -
                             atoms.get_positions()) < 1e-6)
        if format in ['traj', 'cube', 'cfg', 'struct', 'gen']:
            assert np.all(np.abs(a1.get_cell() - atoms.get_cell()) < 1e-6)
        if format in ['cfg']:
            assert np.all(np.abs(a1.get_array('extra') -
                                 atoms.get_array('extra')) < 1e-6)
        if format not in ['cube', 'png', 'eps', 'cfg', 'struct', 'etsf',
                          'gen']:
            a2 = read(fname2)
            a3 = read(fname2, index=0)
            a4 = read(fname2, index=slice(None))
            assert len(a4) == 2
    else:
        print
