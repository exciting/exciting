from ase import Atoms
a = Atoms('H2', positions=[(0, 0, 0), (0, 0, 1.1)])

a.pbc[0] = 1
assert a.pbc.any()
assert not a.pbc.all()
a.pbc = 1
assert a.pbc.all()

a.cell = (1, 2, 3)
a.cell *= 2
a.cell[0, 0] = 3
assert not (a.cell.diagonal() - (3, 4, 6)).any()
