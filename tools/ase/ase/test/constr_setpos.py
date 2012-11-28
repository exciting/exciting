import numpy as np

def array_almost_equal(a1, a2, tol=np.finfo(type(1.0)).eps):
    """Replacement for old numpy.testing.utils.array_almost_equal."""
    return (np.abs(a1 - a2) < tol).all()

from ase.structure import molecule
from ase.constraints import FixAtoms

m = molecule('H2')
c = FixAtoms(indices=[atom.index for atom in m])
m.set_constraint(c)

pos1 = m.get_positions()
# shift z-coordinates by 1.
pos = m.get_positions()
pos[:, 2] += 1.

m.set_positions(pos)
# note that set_positions fails silently to set the new positions
# due to the presence of constraints!
assert array_almost_equal(pos1, m.get_positions())

m.positions = pos
# atoms.positions allows one to set the new positions
# even in the presence of constraints!
assert array_almost_equal(pos, m.get_positions())
