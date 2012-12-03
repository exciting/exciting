"""Lattice data - Face Centered Cubic"""

import numpy as np
from math import sqrt
from symmetry import *

#Definition of symmetries
basesymmetries = [np.matrix([[-1, 0, 0],  #Mirror x-axis
                             [0, 1, 0],
                             [0, 0, 1]]),
                  np.matrix([[1, 0, 0],   #Mirror y-axis
                             [0, -1, 0],
                             [0, 0, 1]]),
                  np.matrix([[1, 0, 0],   #Mirror z-axis
                             [0, 1, 0],
                             [0, 0, -1]]),
                  np.matrix([[1, 0, 0],   #Rotation x-axis (4-fold)
                             [0, 0, -1],
                             [0, 1, 0]]),
                  np.matrix([[0, 0, -1],  #Rotation y-axis (4-fold)
                             [0, 1, 0],
                             [1, 0, 0]]),
                  np.matrix([[0, 1, 0],   #Rotation z-axis (4-fold)
                             [-1, 0, 0],
                             [0, 0, 1]]),
                  np.matrix([[0, 0, 1],   #Rotation (111)-axis (3-fold)
                             [1, 0, 0],
                             [0, 1, 0]]),
                  np.matrix([[0, 0, -1],  #Rotation (11-1)-axis (3-fold)
                             [1, 0, 0],
                             [0, -1, 0]]),
                  np.matrix([[0, 0, 1],   #Rotation (1-11)-axis (3-fold)
                             [-1, 0, 0],
                             [0, -1, 0]]),
                  np.matrix([[0, 0, -1],  #Rotation (-111)-axis (3-fold)
                             [-1, 0, 0],
                             [0, 1, 0]])]

symmetries = get_all_symmetries(basesymmetries, 48)

#Definition of used surfaces
surface_names = [(1,0,0), (-1,0,0),
                 (0,1,0), (0,-1,0),
                 (0,0,1), (0,0,-1),
                 (1,1,0), (-1,-1,0),
                 (1,0,1), (-1,0,-1),
                 (0,1,1), (0,-1,-1),
                 (1,-1,0), (-1,1,0),
                 (1,0,-1), (-1,0,1),
                 (0,1,-1), (0,-1,1),
                 (1,1,1), (-1,-1,-1),
                 (-1,1,1), (1,-1,-1),
                 (1,-1,1), (-1,1,-1),
                 (1,1,-1), (-1,-1,1)]

surface_numbers = {}
for i, s in enumerate(surface_names):
    surface_numbers[s] = i

surface_count = len(surface_names)

surface_mapping = {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4,
                   6: 7, 7: 6, 8: 9, 9: 8, 10: 11, 11: 10,
                   12: 13, 13: 12, 14: 15, 15: 14, 16: 17, 17: 16,
                   18: 19, 19: 18, 20: 21, 21: 20, 22: 23, 23: 22,
                   24: 25, 25: 24}

surface_data = ([{'l': 1.0, 'd': 0.5}] * 6 +
                [{'l': 1.5, 'd': 1.0 / 4.0}] * 12 +
                [{'l': 1.0, 'd': 1.0 / 3.0}] * 8)

surface_symmetries = get_surface_symmetries(symmetries, surface_names, surface_numbers)

def surface_fitting(surfaces):
    for i, n1 in enumerate(np.array(surface_names)):
        d1 = surface_data[i]['d']
        a1 = surfaces[i] * d1
        for j, n2 in enumerate(np.array(surface_names)):
            d2 = surface_data[j]['d']
            a2 = surfaces[j] * d2
            nd = np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2))
            if a2 * nd > a1:
                surfaces[j] = int(round(a1 / (nd * d2)))

    return surfaces

def surface_centering(surfaces, basis='100', debug=0):
    if basis == '100':
        #Centering within the basis {[1,0,0], [0,1,0], [0,0,1]}
        dx = (surfaces[0] - surfaces[1]) // 2
        dy = (surfaces[2] - surfaces[3]) // 2
        dz = (surfaces[4] - surfaces[5]) // 2

        if (dx + dy + dz) % 2 == 1:
            dx -= 1

        if debug:
            print '(%i, %i, %i)' % (dx, dy, dz)
    elif basis == '110':
        #Centering within the basis {[1,1,0], [1,0,1], [0,1,1]}
        dl1 = ((surfaces[6] - surfaces[7]) // 4) * 2
        dl2 = ((surfaces[8] - surfaces[9]) // 4) * 2
        dl3 = ((surfaces[10] - surfaces[11]) // 4) * 2

        #Correction for the none orthogonality of the basis
        t1 = (dl1 != 0 and dl2 != 0 and dl3 == 0)
        t2 = (dl1 != 0 and dl2 == 0 and dl3 != 0)
        t3 = (dl1 != 0 and dl2 == 0 and dl3 != 0)

        if t1 or t2 or t3:
            d1 = (3 * dl1) // 2 - dl2 // 2 - dl3 // 2
            d2 = (3 * dl2) // 2 - dl1 // 2 - dl3 // 2
            d3 = (3 * dl3) // 2 - dl1 // 2 - dl2 // 2
        else:
            d1, d2, d3 = 0, 0, 0

        #Converting to '100' basis
        dx = (d1 + d2) // 2
        dy = (d1 + d3) // 2
        dz = (d2 + d3) // 2

        if debug:
            print ('(%i, %i, %i) -> (%i, %i, %i) -> (%i, %i, %i)' %
                   (dl1, dl2, dl3, d1, d2, d3, dx, dy, dz))
    else:
        dx, dy, dz = 0

    s = np.array(surfaces, int)
    ds = np.array([- dx, dx,
                   - dy, dy,
                   - dz, dz,
                   - dx - dy, dx + dy,
                   - dx - dz, dx + dz,
                   - dy - dz, dy + dz,
                   - dx + dy, dx - dy,
                   - dx + dz, dx - dz,
                   - dy + dz, dy - dz,
                   (-dx - dy - dz) // 2, (dx + dy + dz) // 2,
                   (dx - dy - dz) // 2, (-dx + dy + dz) // 2,
                   (-dx + dy - dz) // 2, (dx - dy + dz) // 2,
                   (-dx - dy + dz) // 2, (dx + dy - dz) // 2], int)

    if (s + ds >= 0).all():
        surfaces = s + ds

    return surfaces

#Definition of the neighbor environment
neighbor_names = [(0.5, 0.5, 0), (-0.5, -0.5, 0),
                  (0.5, 0, 0.5), (-0.5, 0, -0.5),
                  (0, 0.5, 0.5), (0, -0.5, -0.5),
                  (0.5, -0.5, 0), (-0.5, 0.5, 0),
                  (0.5, 0, -0.5), (-0.5, 0, 0.5),
                  (0, 0.5, -0.5), (0, -0.5, 0.5)]

neighbor_numbers = {}
for i, n in enumerate(neighbor_names):
    neighbor_numbers[n] = i

neighbor_positions = np.array(neighbor_names, dtype=float)

neighbor_cutoff = 0.8
neighbor_count = 12

neighbor_mapping = {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4,
                    6: 7, 7: 6, 8: 9, 9: 8, 10: 11, 11: 10}

neighbor_symmetries = get_neighbor_symmetries(symmetries,
                                              neighbor_positions,
                                              neighbor_numbers)

#Definition of the atom types that is used based on the neighborlist
basetype_names = [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                  (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                  (0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1),
                  (0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1),
                  (0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                  (0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1),
                  (0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
                  (0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1),
                  (0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1),
                  (0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0),
                  (1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0),
                  (0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0),
                  (1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0),
                  (0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1),
                  (1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1),
                  (0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1),
                  ]

basetype_data = [{'type': 0,
                  'coordination': 0,
                  'name': 'Free atom, Unknown'},
                 {'type': 1,
                  'coordination': 12,
                  'name': 'bulk'},
                 {'type': 2,
                  'coordination': 8,
                  'name': '100 surface'},
                 {'type': 3,
                  'coordination': 7,
                  'name': '110 surface (top), 111-111 edge, 110-111 edge'},
                 {'type': 4,
                  'coordination': 11,
                  'name': '110 surface (bottom)'},
                 {'type': 5,
                  'coordination': 9,
                  'name': '111 surface'},
                 {'type': 6,
                  'coordination': 6,
                  'name': '100-110 edge'},
                 {'type': 7,
                  'coordination': 7,
                  'name': '100-111 edge'},
                 {'type': 8,
                  'coordination': 6,
                  'name': '111-111-100 corner'},
                 {'type': 9,
                  'coordination': 4,
                  'name': '100 surface ad-atom'},
                 {'type': 10,
                  'coordination': 5,
                  'name': '110 surface ad-atom (bottom), A5 site'},
                 {'type': 11,
                  'coordination': 3,
                  'name': '111 surface ad-atom'},
                 {'type': 12,
                  'coordination': 9,
                  'name': '100 surface with ad-atom'},
                 {'type': 13,
                  'coordination': 8,
                  'name': '110 surface with ad-atom'},
                 {'type': 14,
                  'coordination': 10,
                  'name': '111 surface with ad-atom'},
                 {'type': 15,
                  'coordination': 5,
                  'name': 'B5 site'},
                 ]

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
