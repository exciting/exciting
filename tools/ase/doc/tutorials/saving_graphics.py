# creates:  nice.png

import numpy as np

from ase import Atoms
from ase.io import write

atoms = Atoms('Ag', cell=(2.7, 2.7, 2.7), pbc=True) * (18, 8, 8)

# view with ag
#view(atoms)
rotation = '-70x, -20y, -2z' # found using ag menu 'view -> rotate'

#Make colors
from ase.utils import hsv
colors = hsv(atoms.positions[:, 0])

# Textures
tex = ['jmol',] * 288 + ['glass',] * 288+ ['ase3',] * 288 + ['vmd',] * 288


# keywords
kwargs = { # Keywords that exist for eps, png, and pov
'rotation': rotation,
'show_unit_cell': 2,
'colors': colors,
'radii': None,
}

extra_kwargs = { # For povray files only
'display'      : False, # Display while rendering
'pause'        : False, # Pause when done rendering (only if display)
'transparent'  : False, # Transparent background
'canvas_width' : None,  # Width of canvas in pixels
'canvas_height': None,  # Height of canvas in pixels 
'camera_dist'  : 50.,   # Distance from camera to front atom
'image_plane'  : None,  # Distance from front atom to image plane
                        # (focal depth for perspective)
'camera_type'  : 'perspective', # perspective, ultra_wide_angle
'point_lights' : [],             # [[loc1, color1], [loc2, color2],...]
'area_light'   : [(2., 3., 40.) ,# location
                  'White',       # color
                  .7, .7, 3, 3], # width, height, Nlamps_x, Nlamps_y
'background'   : 'White',        # color
'textures'     : tex, # Length of atoms list of texture names
'celllinewidth': 0.05, # Radius of the cylinders representing the cell
}

# Make flat png file
#write('flat.png', atoms, **kwargs)

# Make the color of the glass beads semi-transparent
colors2 = np.zeros((1152, 4))
colors2[:, :3] = colors
colors2[288: 576, 3] = 0.95
kwargs['colors'] = colors2
kwargs.update(extra_kwargs)

# Make the raytraced image
write('nice.pov', atoms, run_povray=True, **kwargs)
