import os
from ase.structure import molecule
from ase.io.bader import attach_charges

fname = 'ACF.dat'
f = open(fname, 'w')
print >> f, """
   #         X           Y           Z        CHARGE     MIN DIST
 ----------------------------------------------------------------
   1      7.0865      8.5038      9.0672      9.0852      1.3250
   2      7.0865      9.9461      7.9403      0.4574      0.3159
   3      7.0865      7.0615      7.9403      0.4574      0.3159
 ----------------------------------------------------------------
  NUMBER OF ELECTRONS:        9.99999
"""
f.close()

atoms = molecule('H2O')
atoms.set_cell([7.5, 9, 9])
atoms.center()

attach_charges(atoms)
attach_charges(atoms, fname)
os.remove(fname)

for atom in atoms:
    print 'Atom', atom.symbol, 'Bader charge', atom.charge 

