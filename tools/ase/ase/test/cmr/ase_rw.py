import os

from ase.test import NotAvailable

try:
    import cmr
except ImportError:
    raise NotAvailable('CMR is required')

from ase.calculators.emt import EMT
from ase.io import read, write
from ase.structure import molecule

cmr_params = {"db_keywords":["O", "ase"], # keyword
              "molecule":"O2"} #field

m1 = molecule('O2')
m1.set_calculator(EMT())
e1 = m1.get_potential_energy()
write("O2.db", m1, cmr_params = cmr_params)

reread = read("O2.db")
e2 = reread.get_potential_energy()
assert abs(e1 - e2) < 1.e-6, str(e1) + ' ' + str(e2)

db_read = cmr.read("O2.db")
assert "O" in db_read["db_keywords"]
assert "ase" in db_read["db_keywords"]
assert db_read["molecule"] == "O2"

# clean
filename = "O2.db"
if os.path.exists(filename): os.unlink(filename)

