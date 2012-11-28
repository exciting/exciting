from ase import Atom, Atoms
from gpaw import GPAW
a = 4.0
b = a / 2**.5
L = 7.0
al = Atoms([Atom('Al')], cell=(b, b, L), pbc=True)
calc = GPAW(kpts=(4, 4, 1))
al.set_calculator(calc)
al.get_potential_energy()
calc.write('Al100.gpw', 'all')
