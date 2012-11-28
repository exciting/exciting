# creates: Au-wire.png

from ase import Atoms
from ase.io import write

d = 2.9
L = 10.0
wire = Atoms('Au',
             positions=[(0, L / 2, L / 2)],
             cell=(d, L, L),
             pbc=(1, 0, 0))
wire *= (6, 1, 1)
wire.positions[:, 0] -= 2 * d
wire.cell[0, 0] = d
#view(wire, block=1)
write('Au-wire.pov', wire,
      show_unit_cell=2,
      rotation='12x,6y',
      transparent=False, display=False, run_povray=True)
