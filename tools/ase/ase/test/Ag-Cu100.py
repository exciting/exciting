from math import sqrt
from ase import Atom, Atoms
from ase.neb import NEB
from ase.constraints import FixAtoms
from ase.vibrations import Vibrations
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton, BFGS

# Distance between Cu atoms on a (100) surface:
d = 3.6 / sqrt(2)
initial = Atoms('Cu',
                positions=[(0, 0, 0)],
                cell=(d, d, 1.0),
                pbc=(True, True, False))
initial *= (2, 2, 1)  # 2x2 (100) surface-cell

# Approximate height of Ag atom on Cu(100) surfece:
h0 = 2.0
initial += Atom('Ag', (d / 2, d / 2, h0))

if 0:
    view(initial)

# Make band:
images = [initial.copy() for i in range(6)]
neb = NEB(images, climb=True)

# Set constraints and calculator:
constraint = FixAtoms(range(len(initial) - 1))
for image in images:
    image.set_calculator(EMT())
    image.set_constraint(constraint)

# Displace last image:
images[-1].positions[-1] += (d, 0, 0)
#images[-1].positions[-1] += (d, d, 0)

# Relax height of Ag atom for initial and final states:
dyn1 = QuasiNewton(images[0])
dyn1.run(fmax=0.01)
dyn2 = QuasiNewton(images[-1])
dyn2.run(fmax=0.01)

# Interpolate positions between initial and final states:
neb.interpolate()

for image in images:
    print image.positions[-1], image.get_potential_energy()

#dyn = MDMin(neb, dt=0.4)
#dyn = FIRE(neb, dt=0.4)
dyn = BFGS(neb, trajectory='mep.traj')
dyn.run(fmax=0.05)

for image in images:
    print image.positions[-1], image.get_potential_energy()

a = images[0]
vib = Vibrations(a, [4])
vib.run()
print vib.get_frequencies()
vib.summary()
print vib.get_mode(-1)
vib.write_mode(-1, nimages=20)
