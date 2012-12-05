from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize.test import run_test
from gpaw import GPAW

name = 'H2'

def get_atoms():
    cell = (5, 5, 5)
    atoms = Atoms('H2', [(0, 0, 0), (0, 0, 1.4)], cell=cell)
    atoms.center()
    return atoms

def get_calculator_emt():
    calc = EMT()
    return calc

def get_calculator_gpaw():
    calc = GPAW(xc='PBE',txt=None)
    return calc

run_test(get_atoms, get_calculator_emt, name + '-emt')
run_test(get_atoms, get_calculator_gpaw, name + '-gpaw', steps=25)
