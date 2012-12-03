from ase import *
from ase.calculators.abinit import Abinit

a0 = 5.43
bulk = Atoms('Si2', [(0, 0, 0),
                     (0.25, 0.25, 0.25)],
             pbc=True)
b = a0 / 2
bulk.set_cell([(0, b, b),
               (b, 0, b),
               (b, b, 0)], scale_atoms=True)

calc = Abinit(label='Si',
              nbands=8,
              xc='PBE',
              ecut=50 * Ry,
              mix=0.01,
              kpts=[10, 10, 10])

bulk.set_calculator(calc)
e = bulk.get_potential_energy()
