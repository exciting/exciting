# creates: a1.png a2.png a3.png cnt1.png cnt2.png gnr1.png gnr2.png
from ase.io import write
from ase.structure import bulk, nanotube, graphene_nanoribbon
import numpy as np

for i, a in enumerate([
    bulk('Cu', 'fcc', a=3.6),
    bulk('Cu', 'fcc', a=3.6, orthorhombic=True),
    bulk('Cu', 'fcc', a=3.6, cubic=True)]):
    write('a%d.pov' % (i + 1), a,
          show_unit_cell=2, display=False, run_povray=True)

cnt1 = nanotube(6, 0, length=4)
cnt1.rotate('x', 'z', rotate_cell=True)
cnt2 = nanotube(3, 3, length=6, bond=1.4, symbol='Si')
cnt2.rotate('x', 'z', rotate_cell=True)

for i, a in enumerate([cnt1, cnt2]):
    write('cnt%d.pov' % (i + 1), a,
          show_unit_cell=2, display=False, run_povray=True)

ind = [2, 0, 1]
gnr1 = graphene_nanoribbon(3, 4, type='armchair')
gnr1.set_cell(np.diag(gnr1.cell)[ind])
gnr1.positions = gnr1.positions[:, ind]
gnr2 = graphene_nanoribbon(2, 6, type='zigzag', saturated=True,
                           C_H=1.1, C_C=1.4, vacuum=3.0, 
                           magnetic=True, initial_mag=1.12)
gnr2.set_cell(np.diag(gnr2.cell)[ind])
gnr2.positions = gnr2.positions[:, ind]

for i, a in enumerate([gnr1, gnr2]):
    write('gnr%d.pov' % (i + 1), a,
          show_unit_cell=2, display=False, run_povray=True)
