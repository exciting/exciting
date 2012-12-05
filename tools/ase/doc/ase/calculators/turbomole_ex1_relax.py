#!/usr/bin/env python
import os

from ase.io import read,write
from ase.optimize import QuasiNewton

# Turbomole input coordinates must be in the file 'coord'.
# The coordinates are updated to the file 'coord' during the minimization.


test = read('coord')
test.set_calculator(Turbomole())

dyn = QuasiNewton(test, trajectory='test.traj')
dyn.run(fmax=0.01)

write('coord.final.tmol', test)

