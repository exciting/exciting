# creates: transport_setup.png

import numpy as np
from ase import Atoms
from ase.data.molecules import molecule
from ase.io import write

a = 3.92 # Experimental lattice constant
sqrt = np.sqrt
cell = np.array([[a / sqrt(3),          0.,               0.],
                 [         0., a / sqrt(2),               0.],
                 [         0., a / sqrt(8), a * sqrt(3 / 8.)]])
repeat = (1, 3, 3)

A = Atoms('Pt', pbc=True, positions=[(0., 0., 0.)])
B = Atoms('Pt', pbc=True, positions=[(0., 1/3., 1/3.)])
C = Atoms('Pt', pbc=True, positions=[(0., 2/3., 2/3.)])

A *= repeat
B *= repeat
C *= repeat

pyramid_BC = Atoms('Pt4', pbc=True, tags=[1, 1, 1, 2],
                   positions=[( 0., 1/3., 1/3.), # B
                              ( 0., 4/3., 1/3.), # B
                              ( 0., 1/3., 4/3.), # B
                              ( 1., 2/3., 2/3.), # C
                              ])

inv_pyramid_BC = pyramid_BC.copy()
inv_pyramid_BC.positions[:, 0] *= -1

def pos(atoms, x):
    atoms2 = atoms.copy()
    atoms2.translate([x, 0, 0])
    return atoms2

princ = pos(A, 0) + pos(B, 1) + pos(C, 2)
large = (pos(princ, -8) +
         pos(princ, -4) +
         pos(princ, 0) +
         pos(A, 3) +
         pos(pyramid_BC, 4) +
         pos(inv_pyramid_BC, 3) +
         pos(princ, 4) +
         pos(princ, 8))

large.set_cell(cell * repeat, scale_atoms=True)
large.cell[0, 0] = 7 * large.cell[0, 0]

dist=18.
large.cell[0, 0] += dist - cell[0, 0]
large.positions[-(9 * 6 + 4):, 0] += dist - cell[0, 0]

tipL, tipR = large.positions[large.get_tags() == 2]
tipdist = np.linalg.norm(tipL - tipR)

mol = molecule('C6H6', pbc=True, tags=[3] * 6 + [4] * 6)
mol.rotate('y', 'x')
mol.rotate('z', 'y')

large += mol
large.positions[-len(mol):] += tipL
large.positions[-len(mol):, 0] += tipdist / 2

old = large.cell.copy()
large *= (1, 1, 3)
large.set_cell(old)

#view(large)

colors = np.zeros((len(large), 3))
colors[:] = [1., 1., .75]

pr = [.7, .1, .1]
H = [1, 1, 1]
C = [.3, .3, .3]
Pt = [.7, .7, .9]

colors[164:218] = pr # principal layer
colors[289:316] = pr # principal layer
colors[218:289] = Pt # Central region Pt
colors[316:322] = C # Molecule C
colors[322:328] = H # Molecule H


#write('test.png', large, rotation='-90x,-13y', radii=.9, show_unit_cell=0, colors=colors)
write('transport_setup.pov', large, rotation='-90x,-13y', radii=1.06,
      show_unit_cell=0, colors=colors,
      display=False, transparent=False, run_povray=True)
