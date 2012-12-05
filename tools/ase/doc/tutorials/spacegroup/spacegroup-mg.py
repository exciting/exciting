from ase.lattice.spacegroup import crystal

a = 3.21
c = 5.21
mg = crystal('Mg', [(1./3., 2./3., 3./4.)], spacegroup=194,
             cellpar=[a, a, c, 90, 90, 120])
