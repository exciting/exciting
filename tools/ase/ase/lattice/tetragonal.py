"""Function-like objects creating tetragonal lattices.

The following lattice creators are defined:
    SimleTetragonal
    CenteredTetragonal
"""

from ase.lattice.orthorhombic import SimpleOrthorhombicFactory,\
     BodyCenteredOrthorhombicFactory
import numpy as np
from ase.data import reference_states as _refstate


class _Tetragonalize:
    "A mixin class for implementing tetragonal crystals as orthorhombic ones."

    # The name of the crystal structure in ChemicalElements
    xtal_name = "tetragonal"

    def make_crystal_basis(self):
        lattice = self.latticeconstant
        if type(lattice) == type({}):
            lattice['b/a'] = 1.0
        else:
            if len(lattice) == 2:
                lattice = (lattice[0], lattice[0], lattice[1])
            else:
                raise ValueError, "Improper lattice constants for tetragonal crystal."
        self.latticeconstant = lattice
        self.orthobase.make_crystal_basis(self)

class SimpleTetragonalFactory(_Tetragonalize, SimpleOrthorhombicFactory):
    "A factory for creating simple tetragonal lattices."
    orthobase = SimpleOrthorhombicFactory

SimpleTetragonal = SimpleTetragonalFactory()


class CenteredTetragonalFactory(_Tetragonalize,
                                BodyCenteredOrthorhombicFactory):
    "A factory for creating centered tetragonal lattices."
    orthobase = BodyCenteredOrthorhombicFactory

CenteredTetragonal = CenteredTetragonalFactory()
