from math import pi, sqrt
from ase import Atoms, Atom

d = 1.14
a = Atoms([Atom('C', (0, 0, 0)), Atom('O', (d, 0, 0))])
a.rotate_euler(phi=pi/2, theta=pi/4, psi=pi)
for p in a[0].position:
    assert p == 0.0
assert abs(a[1].position[0]) < 1e-15
d2 = d / sqrt(2)
assert abs(a[1].position[1] - d2) < 1e-15
assert abs(a[1].position[2] - d2) < 1e-15
