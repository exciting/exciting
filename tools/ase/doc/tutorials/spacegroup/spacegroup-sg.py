from ase.lattice.spacegroup import Spacegroup

sg = Spacegroup(225)

print 'Space group:', sg.no, sg.symbol
print 'Primitive cell:\n', sg.scaled_primitive_cell
print 'Reciprocal cell:\n', sg.scaled_reciprocal_cell
print 'Lattice type:', sg.lattice
print 'Lattice sub-translations', sg.subtrans
print 'Centerosymmetric', sg.centerosymmetric
print 'Rotation matrices (not including inversions)\n', sg.rotations
print 'Translation vectors (not including inversions)\n', sg.translations
