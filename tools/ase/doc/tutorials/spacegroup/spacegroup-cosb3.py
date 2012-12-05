import numpy as np
import ase.utils.geometry as geometry
import ase.io as io

# Create a new atoms instance with Co at origo including all atoms on the 
# surface of the unit cell
cosb3 = geometry.cut(skutterudite, origo=(0.25, 0.25, 0.25), extend=1.01)

# Define the atomic bonds to show
bondatoms = []
symbols = cosb3.get_chemical_symbols()
for i in xrange(len(cosb3)):
    for j in xrange(i):
        if (symbols[i] == symbols[j] == 'Co' and 
            cosb3.get_distance(i, j) < 4.53):
            bondatoms.append((i, j))
        elif (symbols[i] == symbols[j] == 'Sb' and 
              cosb3.get_distance(i, j) < 2.99):
            bondatoms.append((i, j))

# Create nice-looking image using povray
io.write('spacegroup-cosb3.pov', cosb3,         
         transparent=False, 
         display=False,
         run_povray=True,
         camera_type='perspective',
         canvas_width=320,
         radii=0.4,
         rotation='90y',
         bondlinewidth=0.07,
         bondatoms=bondatoms,
         )
