"""
Function-like objects that creates cubic clusters.
"""

import numpy as np
from ase.cluster.factory import ClusterFactory
from ase.data import reference_states as _refstate

class HexagonalFactory(ClusterFactory):
    spacegroup = 191

    xtal_name = 'hexagonal'

    def get_lattice_constant(self):
        "Get the lattice constant of an element with cubic crystal structure."
        symmetry = _refstate[self.atomic_numbers[0]]['symmetry']
        if symmetry != self.xtal_name:
            raise ValueError, ("Cannot guess the %s " % (self.xtal_name,) +
                               "lattice constant of an element with crystal " +
                               "structure %s." % (symmetry,))
        return _refstate[self.atomic_numbers[0]].copy()

    def set_basis(self):
        lattice = self.lattice_constant
        if isinstance(lattice, dict):
            a = lattice['a']
            try:
                c = lattice['c']
            except KeyError:
                c = a * lattice['c/a']
        else:
            if len(lattice) == 2:
                (a, c) = lattice
            else:
                raise ValueError("Improper lattice constants for %s crystal." % (self.xtal_name,))
        
        self.lattice_constant = (a, c)
        self.lattice_basis = np.array([[a, 0., 0.],
                                       [-a/2., a*np.sqrt(3.)/2., 0.],
                                       [0., 0., c]])
        self.resiproc_basis = self.get_resiproc_basis(self.lattice_basis)

    def set_surfaces_layers(self, surfaces, layers):
        for i, s in enumerate(surfaces):
            if len(s) == 4:
                (a, b, c, d) = s
                if a + b + c != 0:
                    raise ValueError(("(%d,%d,%d,%d) is not a valid hexagonal Miller " +
                                      "index, as the sum of the first three numbers " +
                                      "should be zero.") % (a,b,c,d))
                surfaces[i] = [a, b, d]

        ClusterFactory.set_surfaces_layers(self, surfaces, layers)

Hexagonal = HexagonalFactory()

class HexagonalClosedPackedFactory(HexagonalFactory):
    """A factory for creating HCP clusters."""
    spacegroup = 194

    xtal_name = 'hcp'

    atomic_basis = np.array([[0., 0., 0.],
                             [1./3., 2./3., .5]])

HexagonalClosedPacked = HexagonalClosedPackedFactory()

class GraphiteFactory(HexagonalFactory):
    """A factory for creating graphite clusters."""
    xtal_name = "graphite"

    atomic_basis = np.array([[0., 0., 0.],
                             [1./3., 2./3., 0.],
                             [1./3., 2./3., .5],
                             [2./3., 1./3., .5]])

Graphite = GraphiteFactory()

