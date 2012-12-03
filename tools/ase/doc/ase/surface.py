# creates: fcc100.png fcc110.png bcc100.png fcc111.png bcc110.png bcc111.png hcp0001.png fcc111o.png bcc110o.png bcc111o.png hcp0001o.png ontop-site.png hollow-site.png fcc-site.png hcp-site.png bridge-site.png diamond100.png diamond111.png hcp10m10.png

import os
from ase import Atoms
from ase.io import write
from ase.lattice.surface import (fcc100, fcc110, bcc100, fcc111, 
                                 bcc110, bcc111, hcp0001, hcp10m10,
                                 diamond100, diamond111, add_adsorbate)


surfaces = ['fcc100', 'fcc110', 'bcc100', 'hcp10m10', 'diamond100',
            'fcc111', 'bcc110', 'bcc111', 'hcp0001', 'diamond111']

symbols = {'fcc': 'Cu', 'bcc': 'Fe', 'hcp': 'Ru', 'dia': 'C'}
radii = {'fcc': 1.1, 'bcc': 1.06, 'hcp': 1.08, 'dia': 0.5}
adsorbates= {'ontop': 'H', 'hollow': 'O', 'fcc': 'N', 'hcp': 'C',
             'bridge': 'F'}

def save(name, slab):
    write(name + '.png', slab, show_unit_cell=2, radii=radii[name[:3]],
          scale=10)

for name in surfaces:
    f = eval(name)
    for kwargs in [{}, {'orthogonal': True}]:
        print name, kwargs
        try:
            slab = f(symbols[name[:3]], (3, 4, 5), vacuum=4, **kwargs)
        except (TypeError, NotImplementedError):
            continue
        for site in slab.adsorbate_info['sites']:
            if site.endswith('bridge'):
                h = 1.5
            else:
                h = 1.2
            add_adsorbate(slab, adsorbates.get(site, 'F'), h, site)
        if kwargs:
            name += 'o'
        save(name, slab)

for site, symbol in adsorbates.items():
    write('%s-site.png' % site, Atoms(symbol), radii=1.08, scale=10)
