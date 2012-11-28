"Test that atoms.center() works when adding vacuum ()"

import numpy as np
from math import pi, sqrt, cos
from ase import data
from ase.lattice.cubic import FaceCenteredCubic

def checkang(a, b, phi):
    "Check the angle between two vectors."
    cosphi = np.dot(a,b) / sqrt(np.dot(a,a) * np.dot(b,b))
    assert np.abs(cosphi - cos(phi)) < 1e-10

symb = "Cu"
Z = data.atomic_numbers[symb]
a0 = data.reference_states[Z]['a']

# (100) oriented block
atoms = FaceCenteredCubic(size=(5,5,5), symbol="Cu", pbc=(1,1,0))
assert len(atoms) == 5*5*5*4
c = atoms.get_cell()
checkang(c[0], c[1], pi/2)
checkang(c[0], c[2], pi/2)
checkang(c[1], c[2], pi/2)
assert np.abs(5 * a0 - c[2,2]) < 1e-10

# Add vacuum in one direction
vac = 10.0
atoms.center(axis=2, vacuum=vac)
c = atoms.get_cell()
checkang(c[0], c[1], pi/2)
checkang(c[0], c[2], pi/2)
checkang(c[1], c[2], pi/2)
assert np.abs(4.5 * a0 + 2* vac - c[2,2]) < 1e-10

# Add vacuum in all directions
vac = 4.0
atoms.center(vacuum=vac)
c = atoms.get_cell()
checkang(c[0], c[1], pi/2)
checkang(c[0], c[2], pi/2)
checkang(c[1], c[2], pi/2)
assert np.abs(4.5 * a0 + 2* vac - c[0,0]) < 1e-10
assert np.abs(4.5 * a0 + 2* vac - c[1,1]) < 1e-10
assert np.abs(4.5 * a0 + 2* vac - c[2,2]) < 1e-10

# Now a general unit cell
atoms = FaceCenteredCubic(size=(5,5,5), directions=[[1,0,0], [0,1,0], [1,0,1]],
                          symbol="Cu", pbc=(1,1,0))
assert len(atoms) == 5*5*5*2
c = atoms.get_cell()
checkang(c[0], c[1], pi/2)
checkang(c[0], c[2], pi/4)
checkang(c[1], c[2], pi/2)
assert np.abs(2.5 * a0 - c[2,2]) < 1e-10

# Add vacuum in one direction
vac = 10.0
atoms.center(axis=2, vacuum=vac)
c = atoms.get_cell()
checkang(c[0], c[1], pi/2)
checkang(c[0], c[2], pi/4)
checkang(c[1], c[2], pi/2)
assert np.abs(2 * a0 + 2* vac - c[2,2]) < 1e-10

# Recenter without specifying vacuum
atoms.center()
c = atoms.get_cell()
checkang(c[0], c[1], pi/2)
checkang(c[0], c[2], pi/4)
checkang(c[1], c[2], pi/2)
assert np.abs(2 * a0 + 2* vac - c[2,2]) < 1e-10

# Add vacuum in all directions
vac = 4.0
atoms.center(vacuum=vac)
c = atoms.get_cell()
checkang(c[0], c[1], pi/2)
checkang(c[0], c[2], pi/4)
checkang(c[1], c[2], pi/2)
assert np.abs(4.5 * a0 + 2* vac - c[1,1]) < 1e-10
assert np.abs(2 * a0 + 2* vac - c[2,2]) < 1e-10

