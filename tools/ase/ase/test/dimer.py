from ase import Atom, Atoms
from ase.calculators.lj import LennardJones
from ase.constraints import FixBondLength

dimer = Atoms([Atom('X', (0, 0, 0)),
               Atom('X', (0, 0, 1))],
              calculator=LennardJones(),
              constraint=FixBondLength(0, 1))
print dimer.get_forces()
print dimer.positions
dimer.positions[:] += 0.1
print dimer.positions
dimer.positions[:, 2] += 5.1
print dimer.positions
dimer.positions[:] = [(1,2,3),(4,5,6)]
print dimer.positions
dimer.set_positions([(1,2,3),(4,5,6.2)])
print dimer.positions

