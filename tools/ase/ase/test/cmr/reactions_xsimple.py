import os
import warnings
# cmr calls all available methods in ase.atoms detected by the module inspect.
# Therefore also deprecated methods are called - and we choose to silence those warnings.
warnings.filterwarnings('ignore', 'ase.atoms.*deprecated',)

from ase.test import NotAvailable

# make sure a settings file exist (do this only for the tests please)
try:
    import cmr
except ImportError:
    raise NotAvailable('CMR is required')

from ase.calculators.emt import EMT

from ase.structure import molecule

from ase.io import write

# project id: must uniquely identify the project!
project_id = 'simple reaction energies'

reaction = [('N2', -1), ('N', 2)]

calculator = EMT()

for (formula, coef) in reaction:
    m = molecule(formula)
    m.set_calculator(calculator)
    m.get_potential_energy()
    cmr_params = {
        "db_keywords": [project_id],
        # add project_id also as a field to support search across projects
        "project_id": project_id,
        "formula": formula,
        "calculator": calculator.get_name(),
        }
    write(filename=('reactions_xsimple.%s.db' % formula),
          images=m, format='db', cmr_params=cmr_params)

# analyse the results with CMR

from cmr.ui import DirectoryReader
reader = DirectoryReader('.')

# read all compounds in the project calculated with EMT
all = reader.find(name_value_list=[('calculator', 'EMT')],
                  keyword_list=[project_id])

all.print_table(0, columns=["formula", "ase_potential_energy"])
print

group = cmr.create_group()
group_vars = {"reaction":reaction, "output":"group.db"}
sum = 0.0
for (formula, coef) in reaction:
        data = all.get("formula", formula)
        if data is None:
           print "%s is missing"%formula
           sum = None
           break
        sum += coef*data["ase_potential_energy"]
        group.add(data["db_hash"])

group_vars["result"] = sum
group.write(group_vars)
print "Energy: ",sum
group.dump()

# clean
for (formula, coef) in reaction:
    filename=('reactions_xsimple.%s.db' % formula)
    if os.path.exists(filename): os.unlink(filename)
filename = "group.db"
if os.path.exists(filename): os.unlink(filename)
