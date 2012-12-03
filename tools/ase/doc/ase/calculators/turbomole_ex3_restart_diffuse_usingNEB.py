#!/usr/bin/env python

import os

from ase.io import read
from ase.neb import NEB
from ase.calculators.turbomole import Turbomole
from ase.optimize import BFGS

initial = read('initial.coord')
final = read('final.coord')
os.system('rm -f coord; cp initial.coord coord')

#restart
configs = read('neb.traj@-5:')

band = NEB(configs, climb=True)

#Set calculators
for config in configs:
    config.set_calculator(Turbomole())

# Optimize:
relax = BFGS(band, trajectory='neb.traj')
relax.run(fmax=0.05)

