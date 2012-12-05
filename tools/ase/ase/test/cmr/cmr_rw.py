import os
import warnings
# cmr calls all available methods in ase.atoms detected by the module inspect.
# Therefore also deprecated methods are called - and we choose to silence those warnings.
warnings.filterwarnings('ignore', 'ase.atoms.*deprecated',)

from ase.test import NotAvailable

try:
    import cmr
except ImportError:
    raise NotAvailable('CMR is required')

from ase.calculators.emt import EMT
from ase.structure import molecule

m1 = molecule('O2')
m1.set_calculator(EMT())
e1 = m1.get_potential_energy()

data = cmr.atoms2cmr(m1)
data.set_user_variable("molecule", "O2")
data.set_user_variable("potential", "EMT")
data.set_user_variable("db_keywords", ["O2", "EMT"])

data.write("O2.db")

reread = cmr.read("O2.db")
e2 = reread["ase_potential_energy"]
assert abs(e1-e2) < 1.e-6, str(e1) + ' ' + str(e2)

# clean
filename = "O2.db"
if os.path.exists(filename): os.unlink(filename)
