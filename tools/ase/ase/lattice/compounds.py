"""Function-like objects creating lattices with more than one element.

These lattice creators are mainly intended as examples for how to build you
own.  The following crystal structures are defined:

    B1 = NaCl = Rocksalt
    B2 = CsCl
    B3 = ZnS = Zincblende
    L1_2 = AuCu3
    L1_0 = AuCu
    TRI_Fe2O3
    HEX_Fe2O3
    
"""
from ase.lattice.cubic import FaceCenteredCubicFactory, \
    BodyCenteredCubicFactory, DiamondFactory, SimpleCubicFactory
from ase.lattice.tetragonal import SimpleTetragonalFactory
from ase.lattice.triclinic import TriclinicFactory
from ase.lattice.hexagonal import HexagonalFactory

import numpy as np
from ase.data import reference_states as _refstate


# To prevent a layer of element one on one side, and a layer of
# element two on the other side, NaCl is based on SimpleCubic instead
# of on FaceCenteredCubic
class NaClFactory(SimpleCubicFactory):
    "A factory for creating NaCl (B1, Rocksalt) lattices."

    bravais_basis = [[0, 0, 0], [0, 0, 0.5], [0, 0.5, 0], [0, 0.5, 0.5],
                     [0.5, 0, 0], [0.5, 0, 0.5], [0.5, 0.5, 0],
                     [0.5, 0.5, 0.5]]
    element_basis = (0, 1, 1, 0, 1, 0, 0, 1)
    

B1 = NaCl = Rocksalt = NaClFactory()

class CsClFactory(SimpleCubicFactory):
    "A factory for creating CsCl (B2) lattices."
    bravais_basis = [[0, 0, 0], [0.5, 0.5, 0.5]]
    element_basis = (0, 1)

B2 = CsCl = CsClFactory()


#The zincblende structure is easily derived from Diamond, which
#already has the right basis.
class ZnSFactory(DiamondFactory):
    "A factory for creating ZnS (B3, Zincblende) lattices."
    element_basis = (0, 1)

B3 = ZnS = Zincblende = ZnSFactory()


# The L1_0 structure is "based on FCC", but is a tetragonal distortion
# of fcc.  It must therefore be derived from the base-centered
# tetragonal structure.  That structure, however, does not exist,
# since it is equivalent to a simple tetragonal structure rotated 45
# degrees along the z-axis.  Basing L1_2 on that would however give
# unexpected miller indices.  L1_2 will therefore be based on a simple
# tetragonal structure, but with a basis corresponding to a
# base-centered tetragonal.
class AuCuFactory(SimpleTetragonalFactory):
    "A factory for creating AuCu (L1_0) lattices (tetragonal symmetry)."
    bravais_basis = [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
    element_basis = (0, 1, 1, 0)

AuCu = L1_0 = AuCuFactory()

# The L1_2 structure is "based on FCC", but is really simple cubic
# with a basis.
class AuCu3Factory(SimpleCubicFactory):
    "A factory for creating AuCu3 (L1_2) lattices."
    bravais_basis = [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
    element_basis = (0, 1, 1, 1)

AuCu3 = L1_2 = AuCu3Factory()

class TriclinicFe2O3Factory(TriclinicFactory):
    """A factory for creating hematite (Fe2O3) lattices.

     Rhombohedral unit cell.
     Pauling L, Hendricks S B
     Journal of the American Chemical Society 47 (1925) 781-790

     Example::

         #!/usr/bin/env python
    
         from ase.lattice.hexagonal import *
         from ase.lattice.compounds import *
         import ase.io as io
         from ase import Atoms, Atom
    
         index1=3
         index2=3
         index3=3
         mya = 5.42
         myb = 5.42
         myc = 5.42
         myalpha = 55.28
         mybeta = 55.28
         mygamma = 55.28
         gra = TRI_Fe2O3(symbol = ('Fe', 'O'),
                         latticeconstant={'a':mya,'b':myb, 'c':myc, 
                                          'alpha':myalpha, 
                                          'beta':mybeta,
                                          'gamma':mygamma},
                         size=(index1,index2,index3))
         io.write('rhombohedralUC_Fe2O3.xyz', gra, format='xyz')

     """

    bravais_basis = [[0.10534, 0.10534, 0.10534], [0.39466, 0.39466, 0.39466],
                     [0.60534, 0.60534, 0.60534], [0.89466, 0.89466, 0.89466],
                     [0.30569, 0.69431, 0.00000], [0.69431, 0.00000, 0.30569],
                     [0.00000, 0.30569, 0.69431], [0.19431, 0.80569, 0.50000],
                     [0.80569, 0.50000, 0.19431], [0.50000, 0.19431, 0.80569]]
    element_basis = (0, 0, 0, 0, 1, 1, 1, 1, 1, 1)

TRI_Fe2O3 = TriclinicFe2O3Factory()

class HexagonalFe2O3Factory(HexagonalFactory):
    """A factory for creating hematite (Fe2O3) lattices.
     With hexagonal unit cell.
     Blake R L, Hessevick R E, Zoltai T, Finger L W
     American Mineralogist 51 (1966) 123-129
     5.038 5.038 13.772 90 90 120 R-3c
     Fe       0 0 .3553  .0080  .0080 .00029  .0040      0      0
     O    .3059 0   1/4  .0068  .0083 .00046  .0042 .00058  .0012

     Example:
     #!/usr/bin/env python
     from ase.lattice.hexagonal import *
     from ase.lattice.compounds import *
     import ase.io as io
     from ase import Atoms, Atom

     index1=1
     index2=1
     index3=1
     mya = 5.038
     myb = 5.038
     myc = 13.772
     myalpha = 90
     mybeta = 90
     mygamma = 120
     gra = HEX_Fe2O3(symbol = ('Fe', 'O'),
     latticeconstant={'a':mya,'b':myb, 'c':myc, 
     'alpha':myalpha, 
     'beta':mybeta,
     'gamma':mygamma},
     size=(index1,index2,index3))
     io.write('hexaFe2O3.xyz', gra, format='xyz')
     
     """

    bravais_basis =  [[0.000000, 0.000000, 0.355300],     
                      [0.000000, 0.000000, 0.144700], 
                      [0.000000, 0.000000, 0.644700], 
                      [0.000000, 0.000000, 0.855300],     
                      [0.666667, 0.333333, 0.688633], 
                      [0.666667, 0.333333, 0.478033], 
                      [0.666667, 0.333333, 0.978033],     
                      [0.666667, 0.333333, 0.188633], 
                      [0.333333, 0.666667, 0.021967],  
                      [0.333333, 0.666667, 0.811367],     
                      [0.333333, 0.666667, 0.311367], 
                      [0.333333, 0.666667, 0.521967],
                      #Fe to O here
                      [0.305900, 0.000000, 0.250000],
                      [0.000000, 0.305900, 0.250000],
                      [0.694100, 0.694100, 0.250000],
                      [0.694100, 0.000000, 0.750000],
                      [0.000000, 0.694100, 0.750000],
                      [0.305900, 0.305900, 0.750000],
                      [0.972567, 0.333333, 0.583333],
                      [0.666667, 0.639233, 0.583333],
                      [0.360767, 0.027433, 0.583333],
                      [0.360767, 0.333333, 0.083333],
                      [0.666667, 0.027433, 0.083333],
                      [0.972567, 0.639233, 0.083333],
                      [0.639233, 0.666667, 0.916667],
                      [0.333333, 0.972567, 0.916667],
                      [0.027433, 0.360767, 0.916667],
                      [0.027433, 0.666667, 0.416667],
                      [0.333333, 0.360767, 0.416667],
                      [0.639233, 0.972567, 0.416667]]
    element_basis = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

HEX_Fe2O3 = HexagonalFe2O3Factory()

