from math import sqrt
from ase import Atoms, Atom
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import BFGS, QuasiNewton
from ase.neb import NEB

# Distance between Cu atoms on a (111) surface:
a = 3.6
d = a / sqrt(2)
fcc111 = Atoms(symbols='Cu',
               cell=[(d, 0, 0),
                     (d / 2, d * sqrt(3) / 2, 0),
                     (d / 2, d * sqrt(3) / 6, -a / sqrt(3))],
               pbc=True)
slab = fcc111 * (2, 2, 4)
slab.set_cell([2 * d, d * sqrt(3), 1])
slab.set_pbc((1, 1, 0))
slab.calc = EMT()
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
h = 1.5
slab += Atom('C', (d / 2, -b / 2, h))
slab += Atom('O', (d / 2, +b / 2, h))
s = slab.copy()
dyn = QuasiNewton(slab)
dyn.run(fmax=0.05)
#view(slab)

# Make band:
images = [slab]
for i in range(6):
    image = slab.copy()
    image.set_constraint(constraint)
    image.calc = EMT()
    images.append(image)
image[-2].position = image[-1].position
image[-1].x = d
image[-1].y = d / sqrt(3)
dyn = QuasiNewton(images[-1])
dyn.run(fmax=0.05)
neb = NEB(images, climb=not True)

# Set constraints and calculator:

# Displace last image:

# Relax height of Ag atom for initial and final states:

# Interpolate positions between initial and final states:
neb.interpolate()

for image in images:
    print image.positions[-1], image.get_potential_energy()

#dyn = MDMin(neb, dt=0.4)
#dyn = FIRE(neb, dt=0.01)
dyn = BFGS(neb, maxstep=0.04, trajectory='mep.traj')
#from ase.optimize.oldqn import GoodOldQuasiNewton
#dyn = GoodOldQuasiNewton(neb)
dyn.run(fmax=0.05)

for image in images:
    print image.positions[-1], image.get_potential_energy()

if locals().get('display'):
    import os
    error = os.system('ag mep.traj@-7:')
    assert error == 0
