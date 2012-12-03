"""Function-like object creating hexagonal lattices.

The following lattice creators are defined:

* Hexagonal
* HexagonalClosedPacked
* Graphite
* Graphene

Example for using Graphene to create atoms object gra::

    from ase.lattice.hexagonal import *
    import ase.io as io
    from ase import Atoms, Atom
    
    index1=6
    index2=7
    mya = 2.45
    myc = 20.0
    
    gra = Graphene(symbol = 'C',latticeconstant={'a':mya,'c':myc},
                   size=(index1,index2,1))
    io.write('test.xyz', gra, format='xyz')
"""

from ase.lattice.triclinic import TriclinicFactory
import numpy as np
from ase.data import reference_states as _refstate


class HexagonalFactory(TriclinicFactory):
    "A factory for creating simple hexagonal lattices."
    # The name of the crystal structure in ChemicalElements
    xtal_name = "hexagonal"

    def make_crystal_basis(self):
        "Make the basis matrix for the crystal unit cell and the system unit cell."
        # First convert the basis specification to a triclinic one
        if type(self.latticeconstant) == type({}):
            self.latticeconstant['alpha'] = 90
            self.latticeconstant['beta'] = 90
            self.latticeconstant['gamma'] = 120
            self.latticeconstant['b/a'] = 1.0
        else:
            if len(self.latticeconstant) == 2:
                a,c = self.latticeconstant
                self.latticeconstant = (a,a,c,90,90,120)
            else:
                raise ValueError, "Improper lattice constants for hexagonal crystal."
        TriclinicFactory.make_crystal_basis(self)

    def find_directions(self, directions, miller):
        """Find missing directions and miller indices from the specified ones.

        Also handles the conversion of hexagonal-style 4-index notation to
        the normal 3-index notation.
        """
        directions = list(directions)
        miller = list(miller)
        for obj in (directions,miller):
            for i in range(3):
                if obj[i] is not None:
                    (a,b,c,d) = obj[i]
                    if a + b + c != 0:
                        raise ValueError(
                            ("(%d,%d,%d,%d) is not a valid hexagonal Miller " +
                             "index, as the sum of the first three numbers " +
                             "should be zero.") % (a,b,c,d))
                    x = 4*a + 2*b
                    y = 2*a + 4*b
                    z = 3*d
                    obj[i] = (x,y,z)
        TriclinicFactory.find_directions(self, directions, miller)
        
    def print_directions_and_miller(self, txt=""):
        "Print direction vectors and Miller indices."
        print "Direction vectors of unit cell%s:" % (txt,)
        for i in (0,1,2):
            self.print_four_vector("[]", self.directions[i])
        print "Miller indices of surfaces%s:" % (txt,)
        for i in (0,1,2):
            self.print_four_vector("()", self.miller[i])

    def print_four_vector(self, bracket, numbers):
        bra, ket = bracket
        (x,y,z) = numbers
        a = 2*x - y
        b = -x + 2*y
        c = -x -y
        d = 2*z
        print "   %s%d, %d, %d%s  ~  %s%d, %d, %d, %d%s" % \
              (bra,x,y,z,ket, bra,a,b,c,d,ket)
        
        
Hexagonal = HexagonalFactory()

class HexagonalClosedPackedFactory(HexagonalFactory):
    "A factory for creating HCP lattices."
    xtal_name = "hcp"
    bravais_basis = [[0,0,0], [1.0/3.0, 2.0/3.0, 0.5]]

HexagonalClosedPacked = HexagonalClosedPackedFactory()

class GraphiteFactory(HexagonalFactory):
    "A factory for creating graphite lattices."
    xtal_name = "graphite"
    bravais_basis = [[0,0,0], [1.0/3.0, 2.0/3.0, 0], [1.0/3.0,2.0/3.0,0.5], [2.0/3.0,1.0/3.0,0.5]]

Graphite = GraphiteFactory()

class GrapheneFactory(HexagonalFactory):
    "A factory for creating graphene lattices."
    xtal_name = "graphene"
    bravais_basis = [[0,0,0], [1.0/3.0, 2.0/3.0, 0]]

Graphene = GrapheneFactory()

