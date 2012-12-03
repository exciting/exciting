#!/usr/bin/env python

import os

from ase import Atoms
from ase.calculators.dftb import Dftb
from ase.optimize import QuasiNewton
from ase.io import write

d = 1.10
test = Atoms('2O', positions=[(0., 0., 0.), (0., 0., d)])
test.set_pbc([False, False, False])
test.set_calculator(Dftb(label='o2',write_dftb=True,do_spin_polarized=True,unpaired_electrons=2.0,fermi_temperature=100.0,scc=True))

dyn = QuasiNewton(test, trajectory='test.traj')
dyn.run(fmax=0.01)

write('test.final.xyz', test)

