from ase import Atoms
from ase.dft import STM
from gpaw import GPAW
calc = GPAW('Al100.gpw')
a0 = calc.get_atoms()
stm = STM(calc, [0, 1, 2])
c = stm.get_averaged_current(2.5)
h = stm.scan(c)
print h[8]-h[:, 8]


