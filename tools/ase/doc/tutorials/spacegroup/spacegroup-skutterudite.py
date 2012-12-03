from ase.lattice.spacegroup import crystal

a = 9.04
skutterudite = crystal(('Co', 'Sb'),
                       basis=[(0.25,0.25,0.25), (0.0, 0.335, 0.158)],
                       spacegroup=204, 
                       cellpar=[a, a, a, 90, 90, 90])

#ase.view(skutterudite)
