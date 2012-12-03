"""Function-like object creating triclinic lattices.

The following lattice creator is defined:
    Triclinic
"""

from ase.lattice.bravais import Bravais
import numpy as np
from ase.data import reference_states as _refstate

class TriclinicFactory(Bravais):
    "A factory for creating triclinic lattices."

    # The name of the crystal structure in ChemicalElements
    xtal_name = "triclinic"

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
        "Get the lattice constant of an element with triclinic crystal structure."
        if _refstate[self.atomicnumber]['symmetry'] != self.xtal_name:
            raise ValueError(('Cannot guess the %s lattice constant of'
                              + ' an element with crystal structure %s.')
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
            alpha = lattice['alpha']
            beta = lattice['beta']
            gamma = lattice['gamma']
        else:
            if len(lattice) == 6:
                (a,b,c,alpha,beta,gamma) = lattice
            else:
                raise ValueError, "Improper lattice constants for triclinic crystal."

        degree = np.pi / 180.0
        cosa = np.cos(alpha*degree)
        cosb = np.cos(beta*degree)
        sinb = np.sin(beta*degree)
        cosg = np.cos(gamma*degree)
        sing = np.sin(gamma*degree)
        lattice = np.array([[a,0,0],
                            [b*cosg, b*sing,0],
                            [c*cosb, c*(cosa-cosb*cosg)/sing,
                             c*np.sqrt(sinb**2 - ((cosa-cosb*cosg)/sing)**2)]])
        self.latticeconstant = lattice
        self.miller_basis = lattice
        self.crystal_basis = (self.basis_factor *
                              np.dot(self.int_basis, lattice))
        self.basis = np.dot(self.directions, self.crystal_basis)
        assert abs(np.dot(lattice[0],lattice[1]) - a*b*cosg) < 1e-5
        assert abs(np.dot(lattice[0],lattice[2]) - a*c*cosb) < 1e-5
        assert abs(np.dot(lattice[1],lattice[2]) - b*c*cosa) < 1e-5
        assert abs(np.dot(lattice[0],lattice[0]) - a*a) < 1e-5
        assert abs(np.dot(lattice[1],lattice[1]) - b*b) < 1e-5
        assert abs(np.dot(lattice[2],lattice[2]) - c*c) < 1e-5

Triclinic = TriclinicFactory()
