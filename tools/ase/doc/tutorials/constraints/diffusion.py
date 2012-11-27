# creates:  diffusion-path.png
import os
from ase.io import read, write
if 1:
    execfile('diffusion4.py')
images = [read('mep%d.traj' % i) for i in range(5)]
a = images[0] + images[1] + images[2] + images[3] + images[4]
del a.constraints
a *= (2, 1, 1)
a.set_cell(images[0].get_cell())
write('diffusion-path.pov', a, show_unit_cell=2, rotation='-90x',
      transparent=False, display=False, run_povray=True)
