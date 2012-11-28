from ase.test import NotAvailable
import warnings
# cmr calls all available methods in ase.atoms detected by the module inspect.
# Therefore also deprecated methods are called - and we choose to silence those warnings.
warnings.filterwarnings('ignore', 'ase.atoms.*deprecated',)

try:
    import cmr
except ImportError:
    raise NotAvailable('CMR is required')

from ase.calculators.emt import EMT

from ase.structure import molecule

from ase.io import write
from ase.optimize import QuasiNewton

# see the module for the required format of reactions definition
from ase.test.cmr.reactions import reactions

# assure that all reactions define a reaction_id
for r in reactions:
    assert r[-1][0] == 'reaction_id'

optimize = True

calculator = EMT()

# find names of compounds
# (in one of the most obscure ways - python list flattening)
compounds = [c[0] for c in sum([r[:-1] for r in reactions], [])]
# unique
compounds = list(set(compounds))

for formula in compounds:
    m = molecule(formula)
    m.set_calculator(calculator)
    if optimize:
        dyn = QuasiNewton(m,
                          logfile=('%s.log' % formula),
                          trajectory=('%s.traj' % formula),
                          )
        dyn.run()
    else:
        e = m.get_potential_energy()
        write(filename=('%s.traj' % formula), images=m, format='traj')
