"""Structure data - Hexagonal Closed Packed"""

import numpy as np
from math import sqrt
from symmetry import *

#Definition of symmetries (in hcp basis)
basesymmetries = [np.matrix([[1, 0, 0],   #Mirror y-axis
                             [0, 1, 0],
                             [0, 0, -1]]),
                  np.matrix([[0, 1, 0],   #Rotation z-axix (3-fold)
                             [-1, -1, 0],
                             [0, 0, 1]]),
                  np.matrix([[1, 0, 0],   #Rotation a-axis (2-fold)
                             [-1, -1, 0],
                             [0, 0, -1]]),
                  np.matrix([[-1, -1, 0], #Rotation b-axis (2-fold)
                             [0, 1, 0],
                             [0, 0, -1]]),
                  np.matrix([[0, 1, 0],   #Rotation ab-axis (2-fold)
                             [1, 0, 0],
                             [0, 0, -1]]),
                 ]

symmetries = get_all_symmetries(basesymmetries, 12)

#Definition of surfaces
surface_names = [(0,0,1), (0,0,-1), #(001)
                 (1,-1,0), (-1,1,0), #
                 (2,1,0), (-2,-1,0),
                 (1,2,0), (-1,-2,0),
                 (5,-5,3), (5,-5,-3), #
                 (-5,5,3), (-5,5,-3),
                 (5,10,3), (5,10,-3),
                 (10,5,3), (10,5,-3),
                 (-5,-10,3), (-5,-10,-3),
                 (-10,-5,3), (-10,-5,-3),
                 ]

surface_numbers = {}
for i, s in enumerate(surface_names):
    surface_numbers[s] = i

surface_count = len(surface_names)

surface_mapping = {0:1, 1:0, 2:3, 3:2, 4:5, 5:4, 6:7, 7:6,
                   }

surface_data = ([{'d': 0.5}] * 2 +
                [{'d': 0.5}] * 6 +
                [{'d': 1.0/13.0}] * 12)

surface_symmetries = get_surface_symmetries(symmetries, surface_names, surface_numbers)

#Definition of neighbor environment
neighbor_names = [(1,0,0), (-1,0,0),
                  (0,1,0), (0,-1,0),
                  (1,1,0), (-1,-1,0),
                  (1.0/3.0,2.0/3.0,1.0/2.0), (1.0/3.0,-1.0/3.0,1.0/2.0), (-2.0/3.0,-1.0/3.0,1.0/2.0),
                  (1.0/3.0,2.0/3.0,-1.0/2.0), (1.0/3.0,-1.0/3.0,-1.0/2.0), (-2.0/3.0,-1.0/3.0,-1.0/2.0),
                  (2.0/3.0,1.0/3.0,1.0/2.0), (-1.0/3.0,-2.0/3.0,1.0/2.0), (-1.0/3.0,1.0/3.0,1.0/2.0),
                  (2.0/3.0,1.0/3.0,-1.0/2.0), (-1.0/3.0,-2.0/3.0,-1.0/2.0), (-1.0/3.0,1.0/3.0,-1.0/2.0),
                  ] 

neighbor_numbers = {}
for i, n in enumerate(neighbor_names):
    neighbor_numbers[n] = i

neighbor_positions = np.array(neighbor_names, dtype=float)

neighbor_cutoff = 1.2
neighbor_count = 16

neighbor_mapping = {0:1, 1:0, 2:3, 3:2, 4:5, 5:4, 
                    6:16, 7:17, 8:15, 9:13, 10:14, 11:12,
                    12:11, 13:9, 14:10, 15:8, 16:6, 17:7,
                    }

neighbor_symmetries = get_neighbor_symmetries(symmetries,
                                              neighbor_positions,
                                              neighbor_numbers)

#Definition of the atom types that is used based on the neighborlist
basetype_names = []

basetype_data = []

type_count = len(basetype_names)

type_names = []
type_data = []
for i, n in enumerate(basetype_names):
    type_names.append(n)
    type_data.append(basetype_data[i])

    for sym in neighbor_symmetries:
        new_type = apply_neighbor_symmetry(n, sym)

        if not new_type in type_names:
            type_names.append(new_type)
            type_data.append(basetype_data[i])

type_numbers = {}
for i, n in enumerate(type_names):
    type_numbers[n] = i

#Collect all data
data = {'symmetries': symmetries,
        'surface_names': surface_names,
        'surface_numbers': surface_numbers,
        'surface_data': surface_data,
        'surface_count': surface_count,
        'surface_mapping': surface_mapping,
        'surface_symmetries': surface_symmetries,
        'neighbor_positions': neighbor_positions,
        'neighbor_numbers': neighbor_numbers,
        'neighbor_count': neighbor_count,
        'neighbor_cutoff': neighbor_cutoff,
        'neighbor_mapping': neighbor_mapping,
        'neighbor_symmetries': neighbor_symmetries,
        'type_names': type_names,
        'type_numbers': type_numbers,
        'type_data': type_data,
        'type_count': type_count,
       }

