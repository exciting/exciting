"""Function-like objects creating orthorhombic lattices.

The following lattice creators are defined:
    SimleOrthorhombic
    BaseCenteredOrthorhombic
    BodyCenteredOrthorhombic
    FaceCenteredOrthorhombic
"""

from ase.lattice.bravais import Bravais
import numpy as np
from ase.data import reference_states as _refstate


class SimpleOrthorhombicFactory(Bravais):
    "A factory for creating simple orthorhombic lattices."

    # The name of the crystal structure in ChemicalElements
    xtal_name = "orthorhombic"

    # The natural basis vectors of the crystal structure
    int_basis = np.array([[1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1]])
    basis_factor = 1.0

    # Converts the natural basis back to the crystallographic basis
    inverse_basis = np.array([[1, 0, 0],
                              [0, 1, 0],
                              [0, 0, 1]])
    inverse_basis_factor = 1.0

    def get_lattice_constant(self):
        "Get the lattice constant of an element with orhtorhombic crystal structure."
        if _refstate[self.atomicnumber]['symmetry'] != self.xtal_name:
            raise ValueError, (("Cannot guess the %s lattice constant of"
                                + " an element with crystal structure %s.")
                               % (self.xtal_name,
                                  _refstate[self.atomicnumber]['symmetry']))
        return _refstate[self.atomicnumber].copy()

    def make_crystal_basis(self):
        "Make the basis matrix for the crystal unit cell and the system unit cell."
        lattice = self.latticeconstant
        if type(lattice) == type({}):
            a = lattice['a']
            try:
                b = lattice['b']
            except KeyError:
                b = a * lattice['b/a']
            try:
                c = lattice['c']
            except KeyError:
                c = a * lattice['c/a']
        else:
            if len(lattice) == 3:
                (a,b,c) = lattice
            else:
                raise ValueError, "Improper lattice constants for orthorhombic crystal."

        lattice = np.array([[a,0,0],[0,b,0],[0,0,c]])
        self.latticeconstant = lattice
        self.miller_basis = lattice
        self.crystal_basis = (self.basis_factor *
                              np.dot(self.int_basis, lattice))
        self.basis = np.dot(self.directions, self.crystal_basis)
        self.check_basis_volume()

    def check_basis_volume(self):
        "Check the volume of the unit cell."
        vol1 = abs(np.linalg.det(self.basis))
        vol2 = self.calc_num_atoms() * np.linalg.det(self.latticeconstant)
        if self.bravais_basis is not None:
            vol2 /= len(self.bravais_basis)
        if abs(vol1-vol2) > 1e-5:
            print "WARNING: Got volume %f, expected %f" % (vol1, vol2)

SimpleOrthorhombic = SimpleOrthorhombicFactory()

class BaseCenteredOrthorhombicFactory(SimpleOrthorhombicFactory):
    "A factory for creating base-centered orthorhombic lattices."

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

    def check_basis_volume(self):
        "Check the volume of the unit cell."
        vol1 = abs(np.linalg.det(self.basis))
        vol2 = self.calc_num_atoms() * np.linalg.det(self.latticeconstant) / 2.0
        if abs(vol1-vol2) > 1e-5:
            print "WARNING: Got volume %f, expected %f" % (vol1, vol2)

BaseCenteredOrthorhombic = BaseCenteredOrthorhombicFactory()

class BodyCenteredOrthorhombicFactory(SimpleOrthorhombicFactory):
    "A factory for creating body-centered orthorhombic lattices."

    int_basis = np.array([[-1, 1, 1],
                          [1, -1, 1],
                          [1, 1, -1]])
    basis_factor = 0.5
    inverse_basis = np.array([[0, 1, 1],
                              [1, 0, 1],
                              [1, 1, 0]])
    inverse_basis_factor = 1.0

    def check_basis_volume(self):
        "Check the volume of the unit cell."
        vol1 = abs(np.linalg.det(self.basis))
        vol2 = self.calc_num_atoms() * np.linalg.det(self.latticeconstant) / 2.0
        if abs(vol1-vol2) > 1e-5:
            print "WARNING: Got volume %f, expected %f" % (vol1, vol2)

BodyCenteredOrthorhombic = BodyCenteredOrthorhombicFactory()
class FaceCenteredOrthorhombicFactory(SimpleOrthorhombicFactory):
    "A factory for creating face-centered orthorhombic lattices."

    int_basis = np.array([[0, 1, 1],
                          [1, 0, 1],
                          [1, 1, 0]])
    basis_factor = 0.5
    inverse_basis = np.array([[-1, 1, 1],
                              [1, -1, 1],
                              [1, 1, -1]])
    inverse_basis_factor = 1.0

    def check_basis_volume(self):
        "Check the volume of the unit cell."
        vol1 = abs(np.linalg.det(self.basis))
        vol2 = self.calc_num_atoms() * np.linalg.det(self.latticeconstant) / 4.0
        if abs(vol1-vol2) > 1e-5:
            print "WARNING: Got volume %f, expected %f" % (vol1, vol2)

FaceCenteredOrthorhombic = FaceCenteredOrthorhombicFactory()

