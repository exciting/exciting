from ase.calculators.test import make_test_dft_calculation
from ase.dft.stm import STM

atoms = make_test_dft_calculation()
stm = STM(atoms, [0, 1, 2])
c = stm.get_averaged_current(4.5)
h = stm.scan(c)
