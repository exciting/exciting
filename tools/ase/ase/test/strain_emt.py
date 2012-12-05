"""This test checks that the StrainFilter works using the default
built-in EMT calculator."""

from ase.constraints import StrainFilter
from ase.optimize.mdmin import MDMin
from ase.calculators.emt import EMT
from ase.structure import bulk

cu = bulk('Cu', 'fcc', a=3.6)

cu.set_calculator(EMT(fakestress=True))
f = StrainFilter(cu)
opt = MDMin(f, dt=0.01)
opt.run(0.1, steps=2)
