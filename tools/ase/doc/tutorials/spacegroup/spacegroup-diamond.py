from ase.lattice.spacegroup import crystal

a = 3.57
diamond = crystal('C', [(0,0,0)], spacegroup=227, cellpar=[a, a, a, 90, 90, 90])
