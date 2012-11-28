"""Function-like objects creating cubic lattices (SC, FCC, BCC and Diamond).

The following lattice creators are defined:
    SimpleCubic
    FaceCenteredCubic
    BodyCenteredCubic
    Diamond
"""

from ase.lattice.bravais import Bravais
import numpy as np
from ase.data import reference_states as _refstate

class SimpleCubicFactory(Bravais):
    "A factory for creating simple cubic lattices."

    # The name of the crystal structure in ChemicalElements
    xtal_name = "sc"

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

    # For checking the basis volume
    atoms_in_unit_cell = 1
    
    def get_lattice_constant(self):
        "Get the lattice constant of an element with cubic crystal structure."
        if _refstate[self.atomicnumber]['symmetry'] != self.xtal_name:
            raise ValueError, (("Cannot guess the %s lattice constant of"
                                + " an element with crystal structure %s.")
                               % (self.xtal_name,
                                  _refstate[self.atomicnumber]['symmetry']))
        return _refstate[self.atomicnumber]['a']

    def make_crystal_basis(self):
        "Make the basis matrix for the crystal unit cell and the system unit cell."
        self.crystal_basis = (self.latticeconstant * self.basis_factor
                              * self.int_basis)
        self.miller_basis = self.latticeconstant * np.identity(3)
        self.basis = np.dot(self.directions, self.crystal_basis)
        self.check_basis_volume()

    def check_basis_volume(self):
        "Check the volume of the unit cell."
        vol1 = abs(np.linalg.det(self.basis))
        cellsize = self.atoms_in_unit_cell
        if self.bravais_basis is not None:
            cellsize *= len(self.bravais_basis)
        vol2 = (self.calc_num_atoms() * self.latticeconstant**3 / cellsize)
        assert abs(vol1-vol2) < 1e-5

    def find_directions(self, directions, miller):
        "Find missing directions and miller indices from the specified ones."
        directions = list(directions)
        miller = list(miller)
        # Process keyword "orthogonal"
        self.find_ortho(directions)
        self.find_ortho(miller)
        Bravais.find_directions(self, directions, miller)
        
    def find_ortho(self, idx):
        "Replace keyword 'ortho' or 'orthogonal' with a direction."
        for i in range(3):
            if isinstance(idx[i], str) and (idx[i].lower() == "ortho" or
                                            idx[i].lower() == "orthogonal"):
                if self.debug:
                    print "Calculating orthogonal direction", i
                    print  idx[i-2], "X", idx[i-1],  
                idx[i] = reduceindex(cross(idx[i-2], idx[i-1]))
                if self.debug:
                    print "=", idx[i]
                

SimpleCubic = SimpleCubicFactory()

class FaceCenteredCubicFactory(SimpleCubicFactory):
    "A factory for creating face-centered cubic lattices."

    xtal_name = "fcc"
    int_basis = np.array([[0, 1, 1],
                          [1, 0, 1],
                          [1, 1, 0]])
    basis_factor = 0.5
    inverse_basis = np.array([[-1, 1, 1],
                              [1, -1, 1],
                              [1, 1, -1]])
    inverse_basis_factor = 1.0

    atoms_in_unit_cell = 4
    
FaceCenteredCubic = FaceCenteredCubicFactory()

class BodyCenteredCubicFactory(SimpleCubicFactory):
    "A factory for creating body-centered cubic lattices."

    xtal_name = "bcc"
    int_basis = np.array([[-1, 1, 1],
                          [1, -1, 1],
                          [1, 1, -1]])
    basis_factor = 0.5
    inverse_basis = np.array([[0, 1, 1],
                              [1, 0, 1],
                              [1, 1, 0]])
    inverse_basis_factor = 1.0

    atoms_in_unit_cell = 2
    
BodyCenteredCubic = BodyCenteredCubicFactory()

class DiamondFactory(FaceCenteredCubicFactory):
    "A factory for creating diamond lattices."
    xtal_name = "diamond"
    bravais_basis = [[0,0,0], [0.25, 0.25, 0.25]]
    
Diamond = DiamondFactory()

