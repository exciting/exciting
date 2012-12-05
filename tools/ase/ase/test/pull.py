import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.io import PickleTrajectory

Cu = Atoms('Cu', pbc=(1, 0, 0), calculator=EMT())
traj = PickleTrajectory('Cu.traj', 'w')
for a in np.linspace(2.0, 4.0, 20):
    Cu.set_cell([a, 1, 1], scale_atoms=True)
    traj.write(Cu)
