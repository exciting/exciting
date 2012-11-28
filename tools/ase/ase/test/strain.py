from math import sqrt
from ase import Atoms
from ase.constraints import StrainFilter
from ase.optimize.mdmin import MDMin
from ase.io import PickleTrajectory
try:
    from asap3 import EMT
except ImportError:
    pass
else:
    a = 3.6
    b = a / 2
    cu = Atoms('Cu', cell=[(0,b,b),(b,0,b),(b,b,0)], pbc=1) * (6, 6, 6)

    cu.set_calculator(EMT())
    f = StrainFilter(cu, [1, 1, 1, 0, 0, 0])
    opt = MDMin(f, dt=0.01)
    t = PickleTrajectory('Cu.traj', 'w', cu)
    opt.attach(t)
    opt.run(0.001)

# HCP:
    from ase.lattice.surface import hcp0001
    cu = hcp0001('Cu', (1, 1, 2), a=a / sqrt(2))
    cu.cell[1,0] += 0.05
    cu *= (6, 6, 3)

    cu.set_calculator(EMT())
    f = StrainFilter(cu)
    opt = MDMin(f, dt=0.01)
    t = PickleTrajectory('Cu.traj', 'w', cu)
    opt.attach(t)
    opt.run(0.01)

