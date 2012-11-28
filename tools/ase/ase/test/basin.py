import numpy as np
from math import pi, sqrt
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.optimize.basin import BasinHopping
from ase.io import PickleTrajectory, read
from ase.units import kB

N = 7
R = N**(1./3.)
pos = np.random.uniform(-R, R, (N, 3))
s = Atoms('He' + str(N),
          positions = pos)
s.set_calculator(LennardJones())

ftraj = 'lowest.traj'
traj = PickleTrajectory(ftraj, 'w', s)
bh = BasinHopping(s, 
                  temperature=100 * kB, dr=0.5, 
                  optimizer_logfile=None)
bh.attach(traj)
bh.run(10)

Emin, smin = bh.get_minimum()

# recalc energy
smin.set_calculator(LennardJones())
E = smin.get_potential_energy()
assert abs(E - Emin) < 1e-15
traj.close()
smim = read(ftraj)
E = smin.get_potential_energy()
assert abs(E - Emin) < 1e-15

#view(smin)
