#!/usr/bin/env python

# Copyright 2010 CAMd
# (see accompanying license files for details).
 
from optparse import OptionParser

import numpy as np

from ase.lattice.surface import fcc111, hcp0001, bcc110, diamond111, \
    add_adsorbate
from ase.structure import estimate_lattice_constant
from ase.data import reference_states, atomic_numbers, covalent_radii
from ase.io import write
from ase.visualize import view
from ase.atoms import Atoms, string2symbols
from ase.data.molecules import molecule

def build():
    p = OptionParser(usage='%prog  [options] [ads@]surf [output file]',
                     version='%prog 0.1',
                     description='Example ads/surf: fcc-CO@2x2Ru0001')
    p.add_option('-l', '--layers', type='int',
                 default=4,
                 help='Number of layers.')
    p.add_option('-v', '--vacuum', type='float',
                 default=5.0,
                 help='Vacuum.')
    p.add_option('-x', '--crystal-structure',
                 help='Crystal structure.',
                 choices=['sc', 'fcc', 'bcc', 'hcp'])
    p.add_option('-a', '--lattice-constant', type='float',
                 help='Lattice constant in Angstrom.')
    p.add_option('--c-over-a', type='float',
                 help='c/a ratio.')
    p.add_option('--height', type='float',
                 help='Height of adsorbate over surface.')
    p.add_option('--distance', type='float',
                 help='Distance between adsorbate and nearest surface atoms.')
    p.add_option('-M', '--magnetic-moment', type='float', default=0.0,
                 help='Magnetic moment.')
    p.add_option('-G', '--gui', action='store_true',
                 help="Pop up ASE's GUI.")
    p.add_option('-P', '--python', action='store_true',
                 help="Write Python script.")

    opt, args = p.parse_args()

    if not 1 <= len(args) <= 2:
        p.error("incorrect number of arguments")

    if '@' in args[0]:
        ads, surf = args[0].split('@')
    else:
        ads = None
        surf = args[0]

    if surf[0].isdigit():
        i1 = surf.index('x')
        n = int(surf[:i1])
        i2 = i1 + 1
        while surf[i2].isdigit():
            i2 += 1
        m = int(surf[i1 + 1:i2])
        surf = surf[i2:]
    else:
        n = 1
        m = 1

    if surf[-1].isdigit():
        if surf[1].isdigit():
            face = surf[1:]
            surf = surf[0]
        else:
            face = surf[2:]
            surf = surf[:2]
    else:
        face = None

    Z = atomic_numbers[surf]
    state = reference_states[Z]

    if opt.crystal_structure:
        x = opt.crystal_structure
    else:
        x = state['symmetry']
    
    if opt.lattice_constant:
        a = opt.lattice_constant
    else:
        a = estimate_lattice_constant(surf, x, opt.c_over_a)

    script = ['from ase.lattice.surface import ',
              'vac = %r' % opt.vacuum,
              'a = %r' % a]
    
    if x == 'fcc':
        if face is None:
            face = '111'
        slab = fcc111(surf, (n, m, opt.layers), a, opt.vacuum)
        script[0] += 'fcc111'
        script += ['slab = fcc111(%r, (%d, %d, %d), a, vac)' %
                   (surf, n, m, opt.layers)]
        r = a / np.sqrt(2) / 2
    elif x == 'bcc':
        if face is None:
            face = '110'
        slab = bcc110(surf, (n, m, opt.layers), a, opt.vacuum)
        script[0] += 'bcc110'
        script += ['slab = bcc110(%r, (%d, %d, %d), a, vac)' %
                   (surf, n, m, opt.layers)]
        r = a * np.sqrt(3) / 4
    elif x == 'hcp':
        if face is None:
            face = '0001'
        if opt.c_over_a is None:
            c = np.sqrt(8 / 3.0) * a
        else:
            c = opt.c_over_a * a
        slab = hcp0001(surf, (n, m, opt.layers), a, c, opt.vacuum)
        script[0] += 'hcp0001'
        script += ['c = %r * a' % (c / a),
                   'slab = hcp0001(%r, (%d, %d, %d), a, c, vac)' %
                   (surf, n, m, opt.layers)]
        r = a / 2
    elif x == 'diamond':
        if face is None:
            face = '111'
        slab = diamond111(surf, (n, m, opt.layers), a, opt.vacuum)
        script[0] += 'diamond111'
        script += ['slab = diamond111(%r, (%d, %d, %d), a, vac)' %
                   (surf, n, m, opt.layers)]
        r = a * np.sqrt(3) / 8
    else:
        raise NotImplementedError

    magmom = opt.magnetic_moment
    if magmom is None:
        magmom = {'Ni': 0.6, 'Co': 1.2, 'Fe': 2.3}.get(surf, 0.0)
    slab.set_initial_magnetic_moments([magmom] * len(slab))
    if magmom != 0:
        script += ['slab.set_initial_magnetic_moments([%r] * len(slab))' %
                   magmom]
    
    slab.pbc = 1
    script += ['slab.pbc = True']
    
    name = '%dx%d%s%s' % (n, m, surf, face) 

    if ads:
        site = 'ontop'
        if '-' in ads:
            site, ads = ads.split('-')

        name = site + '-' + ads + '@' + name
        symbols = string2symbols(ads)
        nads = len(symbols) 
        if nads == 1:
            script[:0] = ['from ase import Atoms']
            script += ['ads = Atoms(%r)' % ads]
            ads = Atoms(ads)
        else:
            script[:0] = ['from ase.data.molecules import molecule']
            script += ['ads = molecule(%r)' % ads]
            ads = molecule(ads)

        add_adsorbate(slab, ads, 0.0, site)

        d = opt.distance
        if d is None:
            d = r + covalent_radii[ads[0].number] / 2
        
        h = opt.height
        if h is None:
            R = slab.positions
            y = ((R[:-nads] - R[-nads])**2).sum(1).min()**0.5
            h = (d**2 - y**2)**0.5
        else:
            assert opt.distance is None
        
        slab.positions[-nads:, 2] += h

        script[1] += ', add_adsorbate'
        script += ['add_adsorbate(slab, ads, %r, %r)' % (h, site)]
        
    if len(args) == 2:
        write(args[1], slab)
        script[1:1] = ['from ase.io import write']
        script += ['write(%r, slab)' % args[1]]
    elif not opt.gui:
        write(name + '.traj', slab)
        script[1:1] = ['from ase.io import write']
        script += ['write(%r, slab)' % (name + '.traj')]
        
    if opt.gui:
        view(slab)
        script[1:1] = ['from ase.visualize import view']
        script += ['view(slab)']

    if opt.python:
        print('\n'.join(script))


if __name__ == '__main__':
    build()
