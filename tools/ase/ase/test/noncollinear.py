from ase import Atom, Atoms
a = Atoms('H2')

a[0].magmom = 1
m = a.get_initial_magnetic_moments()
assert m.shape == (2,) and (m == [1, 0]).all()

a[1].magmom = -1
m = a.get_initial_magnetic_moments()
assert m.shape == (2,) and (m == [1, -1]).all()
assert a[1].magmom == -1

a.set_initial_magnetic_moments()
a[0].magmom = (0, 1, 0)
m = a.get_initial_magnetic_moments()
assert m.shape == (2, 3) and (m == [(0, 1, 0), (0, 0, 0)]).all()

a[1].magmom = (1, 0, 0)
m = a.get_initial_magnetic_moments()
assert m.shape == (2, 3) and (m == [(0, 1, 0), (1, 0, 0)]).all()
assert (a[1].magmom == (1, 0, 0)).all()
