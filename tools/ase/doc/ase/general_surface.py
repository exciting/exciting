# creates: s1.png s2.png s3.png s4.png general_surface.pdf
from ase.lattice.surface import surface
s1 = surface('Au', (2, 1, 1), 9)
s1.center(vacuum=10, axis=2)

from ase.lattice import bulk
Mobulk = bulk('Mo', 'bcc', a=3.16, cubic=True)
s2 = surface(Mobulk, (3, 2, 1), 9)
s2.center(vacuum=10, axis=2)

a = 4.0
from ase import Atoms
Pt3Rh = Atoms('Pt3Rh',
              scaled_positions=[(0, 0, 0),
                                (0.5, 0.5, 0),
                                (0.5, 0, 0.5),
                                (0, 0.5, 0.5)],
              cell=[a, a, a],
              pbc=True)
s3 = surface(Pt3Rh, (2, 1, 1), 9)
s3.center(vacuum=10, axis=2)

Pt3Rh.set_chemical_symbols('PtRhPt2')
s4 = surface(Pt3Rh , (2, 1, 1), 9)
s4.center(vacuum=10, axis=2)

from ase.io import write
for atoms, name in [(s1, 's1'), (s2, 's2'), (s3, 's3'), (s4, 's4')]:
    write(name + '.pov', atoms,
          rotation='-90x',
          show_unit_cell=2,
          transparent=False,
          display=False,
          run_povray=True)
    
import os

for i in range(2):
    error = os.system('pdflatex -interaction=nonstopmode general_surface ' +
                      '> /dev/null')
    assert error == 0, 'pdflatex failed'
