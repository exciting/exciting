from math import sqrt, pi
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixBondLengths
from ase.optimize import BFGS, QuasiNewton
from ase.neb import SingleCalculatorNEB
from ase.lattice.surface import fcc111, add_adsorbate
from math import sqrt, cos, sin


zpos = cos(134.3/2.0*pi/180.0)*1.197
xpos = sin(134.3/2.0*pi/180.0)*1.19
co2 = Atoms('COO', positions=[(-xpos+1.2,0,-zpos),
                              (-xpos+1.2,-1.1,-zpos),
                              (-xpos+1.2,1.1,-zpos)])

slab = fcc111('Au', size=(2, 2, 4), vacuum=2*5, orthogonal=True)
slab.center()
add_adsorbate(slab,co2,1.5,'bridge')
slab.set_pbc((True,True,False))
d0 = co2.get_distance(-3, -2)
d1 = co2.get_distance(-3, -1)

calc = EMT()
slab.set_calculator(calc)
constraint = FixBondLengths([[-3,-2],[-3,-1]])
slab.set_constraint(constraint)
dyn = BFGS(slab, trajectory='relax.traj')
dyn.run(fmax=0.05)
assert abs(co2.get_distance(-3, -2) - d0) < 1e-14
assert abs(co2.get_distance(-3, -1) - d1) < 1e-14
