from ase.lattice.spacegroup import crystal

a = 2.87
fe = crystal('Fe', [(0,0,0)], spacegroup=229, cellpar=[a, a, a, 90, 90, 90])
