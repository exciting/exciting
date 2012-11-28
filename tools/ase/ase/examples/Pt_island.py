import numpy as np
from math import sqrt
from ase import Atom, Atoms
from ase.optimize import QuasiNewton, FIRE
from ase.constraints import FixAtoms
from ase.neb import NEB
from ase.io import write, PickleTrajectory
from ase.calculators.emt import ASAP

# Distance between Cu atoms on a (100) surface:
d = 2.74
h1 = d * sqrt(3) / 2
h2 = d * sqrt(2.0 / 3)
initial = Atoms(symbols='Pt',
                positions=[(0, 0, 0)],#(1.37,0.79,2.24),(2.74,1.58,4.48),(0,0,6.72),(1.37,0.79,8.96),(2.74,1.58,11.2)],
                cell=([(d,0,0),(d/2,h1,0),(d/2,h1/3,-h2)]),
                pbc=(True, True, True))
initial *= (7, 8, 6)  # 5x5 (100) surface-cell
cell = initial.get_cell()
cell[2] = (0, 0, 22)
initial.set_cell(cell)
#initial.set_pbc((True,True,False))
# Approximate height of Ag atom on Cu(100) surfece:
h0 = 2.2373
initial += Atom('Pt', (10.96, 11.074, h0))
initial += Atom('Pt', (13.7, 11.074, h0))
initial += Atom('Pt', (9.59, 8.701, h0))
initial += Atom('Pt', (12.33, 8.701, h0))
initial += Atom('Pt', (15.07, 8.701, h0))
initial += Atom('Pt', (10.96, 6.328, h0))
initial += Atom('Pt', (13.7, 6.328, h0))

if 0:
    view(initial)

# Make band:
images = [initial.copy() for i in range(7)]
neb = NEB(images)

# Set constraints and calculator:
indices = np.compress(initial.positions[:, 2] < -5.0, range(len(initial)))
constraint = FixAtoms(indices)
for image in images:
    image.set_calculator(ASAP())
    image.constraints.append(constraint)

# Displace last image:
for i in xrange(1,8,1):
    images[-1].positions[-i] += (d/2, -h1/3, 0)

write('initial.traj', images[0])
# Relax height of Ag atom for initial and final states:
for image in [images[0], images[-1]]:
    QuasiNewton(image).run(fmax=0.01)

if 0:
    write('initial.pckl', image[0])
    write('finial.pckl', image[-1])
# Interpolate positions between initial and final states:
neb.interpolate()

for image in images:
    print image.positions[-1], image.get_potential_energy()

traj = PickleTrajectory('mep.traj', 'w')

dyn = FIRE(neb, dt=0.1)
#dyn = MDMin(neb, dt=0.1)
#dyn = QuasiNewton(neb)
dyn.attach(neb.writer(traj))
dyn.run(fmax=0.01,steps=150)

for image in images:
    print image.positions[-1], image.get_potential_energy()
