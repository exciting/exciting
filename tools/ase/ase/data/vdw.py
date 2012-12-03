# encoding: utf-8
""" Van der Waals radii in [A] taken from
http://www.webelements.com/periodicity/van_der_waals_radius/
and the references given there:
1. A. Bondi, J. Phys. Chem., 1964, 68, 441.
2. L. Pauling, The Nature of the Chemical Bond, 
   Cornell University Press, USA, 1945.
3. J.E. Huheey, E.A. Keiter, and R.L. Keiter in Inorganic Chemistry : 
   Principles of Structure and Reactivity, 4th edition, HarperCollins, 
   New York, USA, 1993.W.W. Porterfield in Inorganic chemistry, 
   a unified approach, Addison Wesley Publishing Co., 
   Reading Massachusetts, USA, 1984.
4. A.M. James and M.P. Lord in Macmillan's Chemical and Physical Data, 
   Macmillan, London, UK, 1992.

additional source from: http://de.wikipedia.org/wiki/Van-der-Waals-Radius
5. Manjeera Mantina, Adam C. Chamberlin, Rosendo Valero, 
   Christopher J. Cramer, Donald G. Truhlar: Consistent van der Waals Radii 
   for the Whole Main Group. In: J. Phys. Chem. A. 2009, 113, 5806–5812, 
   doi:10.1021/jp8111556
"""

import numpy as np

vdw_radii = np.array([
 np.nan, # X
 1.20, # H
 1.40, # He [1]
 1.82, # Li [1]
 1.53, # Be [5]
 1.92, # B [5]
 1.70, # C [1]
 1.55, # N [1]
 1.52, # O [1]
 1.47, # F [1]
 1.54, # Ne [1]
 2.27, # Na [1]
 1.73, # Mg [1]
 1.84, # Al [5]
 2.10, # Si [1]
 1.80, # P [1]
 1.80, # S [1]
 1.75, # Cl [1]
 1.88, # Ar [1]
 2.75, # K [1]
 2.31, # Ca [5]
 np.nan, # Sc
 np.nan, # Ti
 np.nan, # V
 np.nan, # Cr
 np.nan, # Mn
 np.nan, # Fe
 np.nan, # Co
 1.63, # Ni [1]
 1.40, # Cu [1]
 1.39, # Zn [1]
 1.87, # Ga [1]
 2.11, # Ge [5]
 1.85, # As [1]
 1.90, # Se [1]
 1.85, # Br [1]
 2.02, # Kr [1]
 3.03, # Rb [5]
 2.49, # Sr [5]
 np.nan, # Y
 np.nan, # Zr
 np.nan, # Nb
 np.nan, # Mo
 np.nan, # Tc
 np.nan, # Ru
 np.nan, # Rh
 1.63, # Pd [1]
 1.72, # Ag [1]
 1.58, # Cd [1]
 1.93, # In [1]
 2.17, # Sn [1]
 2.06, # Sb [5]
 2.06, # Te [1]
 1.98, # I [1]
 2.16, # Xe [1]
 3.43, # Cs [5]
 2.49, # Ba [5]
 np.nan, # La
 np.nan, # Ce
 np.nan, # Pr
 np.nan, # Nd
 np.nan, # Pm
 np.nan, # Sm
 np.nan, # Eu
 np.nan, # Gd
 np.nan, # Tb
 np.nan, # Dy
 np.nan, # Ho
 np.nan, # Er
 np.nan, # Tm
 np.nan, # Yb
 np.nan, # Lu
 np.nan, # Hf
 np.nan, # Ta
 np.nan, # W
 np.nan, # Re
 np.nan, # Os
 np.nan, # Ir
 1.75, # Pt [1]
 1.66, # Au [1]
 1.55, # Hg [1]
 1.96, # Tl [1]
 2.02, # Pb [1]
 2.07, # Bi [5]
 1.97, # Po [5]
 2.02, # At [5]
 2.20, # Rn [5]
 3.48, # Fr [5]
 2.83, # Ra [5]
 np.nan, # Ac
 np.nan, # Th
 np.nan, # Pa
 1.86, # U [1]
 np.nan, # Np
 np.nan, # Pu
 np.nan, # Am
 np.nan, # Cm
 np.nan, # Bk
 np.nan, # Cf
 np.nan, # Es
 np.nan, # Fm
 np.nan, # Md
 np.nan, # No
 np.nan]) # Lr
