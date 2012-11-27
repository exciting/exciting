# -*- coding: utf-8 -*-
# creates:  surface.png
import os
from ase.io import read, write
execfile('N2Cu.py')
image = read('N2Cu.traj@-1')
write('surface.pov', image, transparent=False, display=False, run_povray=True)
