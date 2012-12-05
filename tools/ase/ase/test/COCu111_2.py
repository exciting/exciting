from math import sqrt
from ase import Atoms, Atom
from ase.constraints import FixAtoms
from ase.optimize import FIRE, QuasiNewton, BFGS
from ase.neb import SingleCalculatorNEB
from ase.calculators.emt import EMT

Optimizer = BFGS

# Distance between Cu atoms on a (111) surface:
a = 3.6
d = a / sqrt(2)
fcc111 = Atoms(symbols='Cu',
               cell=[(d, 0, 0),
                     (d / 2, d * sqrt(3) / 2, 0),
                     (d / 2, d * sqrt(3) / 6, -a / sqrt(3))],
               pbc=True)
initial = fcc111 * (2, 2, 4)
initial.set_cell([2 * d, d * sqrt(3), 1])
initial.set_pbc((1, 1, 0))
initial.set_calculator(EMT())
Z = initial.get_positions()[:, 2]
indices = [i for i, z in enumerate(Z) if z < Z.mean()]
constraint = FixAtoms(indices=indices)
initial.set_constraint(constraint)
dyn = Optimizer(initial)
dyn.run(fmax=0.05)
Z = initial.get_positions()[:, 2]
print Z[0] - Z[1]
print Z[1] - Z[2]
print Z[2] - Z[3]

b = 1.2
h = 1.5
initial += Atom('C', (d / 2, -b / 2, h))
initial += Atom('O', (d / 2, +b / 2, h))
s = initial.copy()
dyn = Optimizer(initial)
dyn.run(fmax=0.05)
#view(initial)

# create final
final = initial.copy()
final.set_calculator(EMT())
final.set_constraint(constraint)
final[-2].position = final[-1].position
final[-1].x = d
final[-1].y = d / sqrt(3)
dyn = Optimizer(final)
dyn.run(fmax=0.1)
#view(final)

# create 2 intermediate step neb
neb = SingleCalculatorNEB([initial, final])
neb.refine(2)
neb.set_calculators(EMT())
assert neb.n() == 4

dyn = Optimizer(neb, maxstep=0.04, trajectory='mep_2coarse.traj')
dyn.run(fmax=0.1)
#dyn.run(fmax=39.1)

# read from the trajectory
neb = SingleCalculatorNEB('mep_2coarse.traj@-4:')

# refine in the important region
neb.refine(2, 1, 3)
neb.set_calculators(EMT())
dyn = Optimizer(neb, maxstep=0.04, trajectory='mep_2fine.traj')
dyn.run(fmax=0.1)
assert len(neb.images) == 8
