"""Function-like object creating monoclinic lattices.

The following lattice creator is defined:
    SimpleMonoclinic
    BaseCenteredMonoclinic
"""

from ase.lattice.triclinic import TriclinicFactory
import numpy as np
from ase.data import reference_states as _refstate


class SimpleMonoclinicFactory(TriclinicFactory):
    "A factory for creating simple monoclinic lattices."
    # The name of the crystal structure in ChemicalElements
    xtal_name = "monoclinic"

    def make_crystal_basis(self):
        "Make the basis matrix for the crystal unit cell and the system unit cell."
        # First convert the basis specification to a triclinic one
        if type(self.latticeconstant) == type({}):
            self.latticeconstant['beta'] = 90
            self.latticeconstant['gamma'] = 90
        else:
            if len(self.latticeconstant) == 4:
                self.latticeconstant = self.latticeconstant + (90,90)
            else:
                raise ValueError, "Improper lattice constants for monoclinic crystal."

        TriclinicFactory.make_crystal_basis(self)
        
SimpleMonoclinic = SimpleMonoclinicFactory()

class BaseCenteredMonoclinicFactory(SimpleMonoclinicFactory):
    # The natural basis vectors of the crystal structure
    int_basis = np.array([[1, -1, 0],
                          [1, 1, 0],
                          [0, 0, 2]])
    basis_factor = 0.5

    # Converts the natural basis back to the crystallographic basis
    inverse_basis = np.array([[1, 1, 0],
                              [-1, 1, 0],
                              [0, 0, 1]])
    inverse_basis_factor = 1.0

BaseCenteredMonoclinic = BaseCenteredMonoclinicFactory()
