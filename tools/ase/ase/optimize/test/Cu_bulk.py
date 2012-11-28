from ase.calculators.emt import EMT
from ase.lattice.cubic import FaceCenteredCubic
from ase.optimize.test import run_test

name = 'Cu_bulk'

def get_atoms():
    atoms = FaceCenteredCubic(directions=[[1,-1,0], [1,1,0], [0,0,1]],
                              size=(3,3,3), symbol='Cu', pbc=(1,1,1))
    atoms.rattle(stdev=0.1,seed=42)
    return atoms

def get_calculator():
    return EMT()

run_test(get_atoms, get_calculator, name, fmax=0.02)
