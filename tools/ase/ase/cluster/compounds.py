import numpy as np

from ase.cluster.cubic import SimpleCubicFactory

# The L1_2 structure is "based on FCC", but is really simple cubic
# with a basis.
class AuCu3Factory(SimpleCubicFactory):
    "A factory for creating AuCu3 (L1_2) lattices."

    atomic_basis = np.array([[0., 0., 0.],
                             [0., .5, .5],
                             [.5, 0., .5],
                             [.5, .5, 0.]])

    element_basis = [0, 1, 1, 1]

AuCu3 = L1_2 = AuCu3Factory()

