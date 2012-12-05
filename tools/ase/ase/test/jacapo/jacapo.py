# do some tests here before we import
# Right version of Scientific?
from ase.test import NotAvailable
import os
try:
    import Scientific
    version = Scientific.__version__.split(".")
    print 'Found ScientificPython version: ',Scientific.__version__
    if map(int,version) < [2,8]:
        print 'ScientificPython 2.8 or greater required for numpy support in NetCDF'
        raise NotAvailable('ScientificPython version 2.8 or greater is required')
except (ImportError, NotAvailable):
    print "No Scientific python found. Check your PYTHONPATH"
    raise NotAvailable('ScientificPython version 2.8 or greater is required')

if not (os.system('which dacapo.run') == 0):
    print "No Dacapo Fortran executable (dacapo.run) found. Check your path settings."
    raise NotAvailable('dacapo.run is not installed on this machine or not in the path')

# Now Scientific 2.8 and dacapo.run should both be available

from ase import Atoms, Atom
from ase.calculators.jacapo import Jacapo

atoms = Atoms([Atom('H',[0,0,0])],
            cell=(2,2,2))

calc = Jacapo('Jacapo-test.nc',
              pw=200,
              nbands=2,
              kpts=(1,1,1),
              spinpol=False,
              dipole=False,
              symmetry=False,
              ft=0.01)

atoms.set_calculator(calc)

print atoms.get_potential_energy()
os.system('rm -f Jacapo-test.nc Jacapo-test.txt')
