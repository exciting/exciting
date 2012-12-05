import numpy as np
from ase import Atoms
from ase.units import fs
from ase.calculators.test import TestPotential
from ase.calculators.emt import EMT
from ase.md import VelocityVerlet
from ase.io import PickleTrajectory, read
from ase.optimize import QuasiNewton

np.seterr(all='raise')
a = Atoms('4X',
          masses=[1, 2, 3, 4],
          positions=[(0, 0, 0),
                     (1, 0, 0),
                     (0, 1, 0),
                     (0.1, 0.2, 0.7)],
          calculator=TestPotential())
print a.get_forces()
md = VelocityVerlet(a, dt=0.5 * fs, logfile='-', loginterval=500)
traj = PickleTrajectory('4N.traj', 'w', a)
md.attach(traj.write, 100)
e0 = a.get_total_energy()
md.run(steps=10000)
del traj
assert abs(read('4N.traj').get_total_energy() - e0) < 0.0001

qn = QuasiNewton(a)
qn.run(0.001)
assert abs(a.get_potential_energy() - 1.0) < 0.000002
