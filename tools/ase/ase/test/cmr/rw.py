import os
import warnings
# cmr calls all available methods in ase.atoms detected by the module inspect.
# Therefore also deprecated methods are called - and we choose to silence those warnings.
warnings.filterwarnings('ignore', 'ase.atoms.*deprecated',)

import numpy as np

def array_almost_equal(a1, a2, tol=np.finfo(type(1.0)).eps):
    """Replacement for old numpy.testing.utils.array_almost_equal."""
    return (np.abs(a1 - a2) < tol).all()

from ase.test import NotAvailable

# make sure a settings file exist (do this only for the tests please)
# this test should be run with cmr!
try:
    import cmr
except ImportError:
    raise NotAvailable('CMR is required')

from ase.calculators.emt import EMT

from ase.io import read, write

from ase.structure import molecule

m1 = molecule('O2')
m1.center(2.0)

write("O2.db", images=m1)

m1.set_calculator(EMT())
e1 = m1.get_potential_energy()
f1 = m1.get_forces()

m2 = read("O2.db")

m2.set_calculator(EMT())
e2 = m2.get_potential_energy()
f2 = m1.get_forces()

# assume atoms definitions are the same if energy/forces are the same: can we do better?
assert abs(e1-e2) < 1.e-6, str(e1) + ' ' + str(e2)
assert array_almost_equal(f1, f2, tol=1.e-6)

# clean
filename = "O2.db"
if os.path.exists(filename): os.unlink(filename)

