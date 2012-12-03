"""Test that atoms.get_center_of_mass(scaled=True) works"""

import numpy as np
from ase import Atoms

d = 1.142
a = Atoms('CO', positions=[(2, 0, 0), (2, -d, 0)], pbc=True)
a.set_cell(np.array(((4, -4, 0), (0, 5.657, 0), (0, 0, 10))))

def array_almost_equal(a1, a2, tol=np.finfo(type(1.0)).eps):
    return (np.abs(a1 - a2) < tol).all()

scaledref = np.array((0.5, 0.23823499, 0.))
assert array_almost_equal(a.get_center_of_mass(scaled=True), scaledref, tol=1e-8)
