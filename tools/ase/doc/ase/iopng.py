# creates: io1.png io2.png io3.png

from ase import Atoms
from ase.lattice.surface import fcc111, add_adsorbate
from ase.io import write, read

adsorbate = Atoms('CO')
adsorbate[1].z = 1.1
a = 3.61
slab = fcc111('Cu', (2, 2, 3), a=a, vacuum=7.0)
add_adsorbate(slab, adsorbate, 1.8, 'ontop')

#view(slab)
write('io1.png', slab * (3, 3, 1), rotation='10z,-80x')
write('io2.pov', slab * (3, 3, 1), rotation='10z,-80x',
      transparent=False, display=False, run_povray=True)
d = a / 2**0.5
write('io3.pov', slab * (2, 2, 1),
      bbox=(d, 0, 3 * d, d * 3**0.5),
      transparent=False, display=False, run_povray=True)

write('slab.xyz', slab)
a = read('slab.xyz')
a.get_cell()
a.get_pbc()
write('slab.traj', slab)
b = read('slab.traj')
b.get_cell()
b.get_pbc()
