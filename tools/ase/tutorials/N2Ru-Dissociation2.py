import numpy as np

from ase.io import read
from ase.constraints import FixAtoms
from ase.calculators.emt import EMT
from ase.neb import NEB
from ase.optimize import QuasiNewton

initial = read('N2.traj')
final = read('2N.traj')

configs = [initial.copy() for i in range(8)] + [final]

constraint = FixAtoms(mask=[atom.symbol != 'N' for atom in initial])
for config in configs:
    config.set_calculator(EMT())
    config.set_constraint(constraint)
    
band = NEB(configs)
band.interpolate()

# Create a quickmin object:
relax = QuasiNewton(band)

relax.run(steps=20)

e0 = initial.get_potential_energy()
for config in configs:
    d = config[-2].position - config[-1].position
    print np.linalg.norm(d), config.get_potential_energy() - e0
