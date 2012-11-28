from math import sin, cos, pi
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize.test import run_test
from ase.lattice.surface import fcc111, add_adsorbate

name = 'CO_Au111'

def get_atoms():
    zpos = cos(134.3/2.0*pi/180.0)*1.197
    xpos = sin(134.3/2.0*pi/180.0)*1.19
    no2 =Atoms('CO', positions=[(-xpos+1.2,0,-zpos), (-xpos+1.2,-1.1,-zpos)])

    # Surface slab
    slab =fcc111('Au', size=(2, 2, 4),vacuum=2*5, orthogonal = True )
    slab.center()
    add_adsorbate(slab,no2,1.5,'bridge')
    slab.set_pbc((True,True,False))

    #constraints
    constraint = FixAtoms(mask=[(a.tag == 4) or (a.tag == 3) or (a.tag==2) for a in slab])
    slab.set_constraint(constraint)
    return slab

def get_calculator():
    calc = EMT()
    return calc

run_test(get_atoms, get_calculator, name, steps=200)
