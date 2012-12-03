from ase import Atoms, Atom
from ase.calculators.vasp import Vasp

a = [6.5, 6.5, 7.7]
d = 2.3608
NaCl = Atoms([Atom('Na', [0, 0, 0], magmom=1.928),
              Atom('Cl', [0, 0, d], magmom=0.75)],
             cell=a)

calc = Vasp(prec = 'Accurate', 
            xc = 'PBE', 
            lreal = False)
NaCl.set_calculator(calc)

print NaCl.get_magnetic_moment()
