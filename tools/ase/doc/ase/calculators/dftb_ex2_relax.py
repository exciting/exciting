#!/usr/bin/python

from ase.io import read, write
from ase.calculators.dftb import Dftb
from ase.optimize import QuasiNewton

system = read('h2o_1.xyz')
system.set_calculator(Dftb(label='dftb',write_dftb=False))

qn = QuasiNewton(system)
qn.run(fmax=0.005)

write('final.xyz', system)
write('final.traj', system)
