# creates: a1.png a2.png a3.png

import numpy as np

from ase import Atom,Atoms
from ase.io import write

atoms = Atoms('Ni4', positions= [(0, 0, 0),
                                 (0.45, 0, 0),
                                 (0, 0.5, 0),
                                 (0.5, 0.5, 0),])
atoms[1].x = 0.5
a = 3.55
cell = [(2 / np.sqrt(2.0) * a, 0, 0),
        (1 / np.sqrt(2.0) * a, np.sqrt(3.0 / 2.0) * a, 0),
        (0, 0, 10 * np.sqrt(3.0) / 3.0 * a)]
atoms.set_cell(cell, scale_atoms=True)
write('a1.png', atoms, rotation='-73x', show_unit_cell=2)

a = atoms.repeat((3, 3, 2))
a.set_cell(atoms.get_cell())
write('a2.png', a, rotation='-73x', show_unit_cell=True)

xyzcell = np.identity(3) # The 3x3 unit matrix
atoms.set_cell(xyzcell, scale_atoms=True)  # Set the unit cell and rescale
atoms.append(Atom('Ni', (1/6., 1/6., .1)))  
atoms.set_cell(cell, scale_atoms=True)  # Set the unit cell and scale back
a = atoms.repeat((3, 3, 1))
a.set_cell(atoms.get_cell())
write('a3.png', a, show_unit_cell=True)

