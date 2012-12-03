from math import sqrt
from ase import Atoms, Atom
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.io import PickleTrajectory
from ase.neb import NEB
from ase.calculators.emt import EMT

# Distance between Cu atoms on a (111) surface:
a = 3.6
d = a / sqrt(2)
y = d * sqrt(3) / 2
fcc111 = Atoms('Cu',
               cell=[(d, 0, 0),
                     (d / 2, y, 0),
                     (d / 2, y / 3, -a / sqrt(3))],
               pbc=True)
slab = fcc111 * (2, 2, 4)
slab.set_cell([2 * d, 2 * y, 1])
slab.set_pbc((1, 1, 0))
slab.set_calculator(EMT())
Z = slab.get_positions()[:, 2]
indices = [i for i, z in enumerate(Z) if z < Z.mean()]
constraint = FixAtoms(indices=indices)
slab.set_constraint(constraint)
dyn = QuasiNewton(slab)
dyn.run(fmax=0.05)
Z = slab.get_positions()[:, 2]
print Z[0] - Z[1]
print Z[1] - Z[2]
print Z[2] - Z[3]

b = 1.2
h = 2.0
slab += Atom('C', (d, 2 * y / 3, h))
slab += Atom('O', (3 * d / 2, y / 3, h))
traj = PickleTrajectory('initial.traj', 'w', slab)
dyn = QuasiNewton(slab)
dyn.attach(traj.write)
dyn.run(fmax=0.05)
#view(slab)
# Make band:
images = [slab.copy() for i in range(6)]
neb = NEB(images, climb=True)

# Set constraints and calculator:
for image in images:
    image.set_calculator(EMT())
    image.set_constraint(constraint)

# Displace last image:
images[-1].positions[-1] = (2 * d, 2 * y / 3, h)
traj = PickleTrajectory('final.traj', 'w', images[-1])
dyn = QuasiNewton(images[-1])
dyn.attach(traj.write)
dyn.run(fmax=0.05)

# Interpolate positions between initial and final states:
neb.interpolate()

for image in images:
    print image.positions[-1], image.get_potential_energy()

traj = PickleTrajectory('mep.traj', 'w')

#dyn = MDMin(neb, dt=0.4)
#dyn = FIRE(neb, dt=0.4)
dyn = QuasiNewton(neb)
dyn.attach(neb.writer(traj))
dyn.run(fmax=0.05)

for image in images:
    print image.positions[-1], image.get_potential_energy()
