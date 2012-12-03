from fcc import *
import numpy as np

def get_all_symmetries(symmetries=None, max=99):
    if symmetries is None:
        raise Warning('Some unique symmetries are needed to start.')

    symmetries_all = symmetries[:]

    for i, l in enumerate(symmetries):
        for j, m, in enumerate(symmetries):
            if len(symmetries_all) == max: break
            v = l * m

            exist = False
            for w in symmetries_all:
                if (v == w).all():
                    exist = True
                    break

            if not exist:
                #print 'New added: %i x %i' % (i, j)
                symmetries_all.append(v)

    for i, l in enumerate(symmetries):
        for j, m, in enumerate(symmetries):
            for k, n in enumerate(symmetries):
                if len(symmetries_all) == max: break
                v = l * m * n

                exist = False
                for w in symmetries_all:
                    if (v == w).all():
                        exist = True
                        break

                if not exist:
                    #print 'New added: %i x %i x %i' % (i, j, k)
                    symmetries_all.append(v)

    #print 'There are %i symmetry operations.' % len(symmetries_all)
    return symmetries_all

def get_neighbor_symmetries(symmetries=None, neighbor_positions=None, neighbor_numbers=None):
    if symmetries is None or neighbor_positions is None or neighbor_numbers is None:
        raise Warning('Both symmetries, positions and numbers for the neighbors are needed.')

    neighbor_symmetries = []

    for s in symmetries:
        neighbor_symmetry = []
        
        for p in neighbor_positions:
            new_p = np.array(np.matrix(p) * s)
            neighbor_symmetry.append(neighbor_numbers[tuple(new_p[0])])

        neighbor_symmetries.append(neighbor_symmetry)

    return neighbor_symmetries

def get_surface_symmetries(symmetries=None, surface_names=None, surface_numbers=None):
    if symmetries is None or surface_names is None or surface_numbers is None:
        raise Warning('Both symmetries, names and numbers for the surfaces are needed.')

    surface_symmetries = []

    for sym in symmetries:
        surface_symmetry = []

        for s in surface_names:
            ns = np.array(np.matrix(s) * sym)
            surface_symmetry.append(surface_numbers[tuple(ns[0])])

        surface_symmetries.append(surface_symmetry)

    return np.array(surface_symmetries, int)

def apply_neighbor_symmetry(neighbors=None, symmetry=None):
    if neighbors is None or symmetry is None:
        raise Warning('Both neighbor list and symmetry list are needed.')

    new_neighbors = [0] * len(symmetry)
    
    for i, n in enumerate(symmetry):
        new_neighbors[i] = neighbors[n]

    return tuple(new_neighbors)
