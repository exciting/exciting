from ase import Atoms, Atom, view
from gpaw import GPAW

logo = """\
 H   HH HHH
H H H   H
HHH  H  HH
H H   H H
H H HH  HHH"""

d = 0.8
atoms = Atoms()
for y, line in enumerate(logo.split('\n')):
    for x, c in enumerate(line):
        if c == 'H':
            atoms.append(Atom('H', [d * x, -d * y, 0]))
atoms.center(vacuum=2.0)
view(atoms)

if 0:
    calc = GPAW(nbands=30)
    atoms.set_calculator(calc)
    atoms.get_potential_energy()
    calc.write('ase-logo.gpw')
else:
    calc = GPAW('ase-logo.gpw', txt=None)

density = calc.get_pseudo_density()
image = density[..., density.shape[2] // 2]

if 1: # scale colors to wiki background / foreground
    import numpy as np
    background = np.array([[[19., 63., 82.]]]).T / 255 # 1c4e63 blueish
    foreground = np.array([[[1., 1., 1.]]]).T  # white
    image = background + image / image.max() * (foreground - background)
    image = image.T
else: # Use a standard color scheme
    image = pl.cm.hot(image.T)

import pylab as pl
x, y, z = atoms.cell.diagonal()
pl.figure(1, figsize=(8, 8 * y / x), dpi=80)
pl.axes([0, 0, 1, 1])
pl.imshow(image,
          origin='lower',
          extent=[0, x, 0, y],
          aspect='equal',
          interpolation='spline16')
pl.axis('off')
pl.savefig('ase-logo.png', dpi=80)
pl.show()
