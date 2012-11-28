import numpy as np

from ase import Atoms
from ase.data import atomic_numbers, reference_states

def Icosahedron(symbol, noshells, latticeconstant=None):
    """
    Returns a cluster with the icosahedra symmetry.

    Parameters
    ----------
    symbol: The chemical symbol (or atomic number) of the element.

    noshells: The number of shells (>= 1).

    latticeconstant (optional): The lattice constant. If not given,
    then it is extracted form ase.data.
    """

    # Interpret symbol
    if isinstance(symbol, str):
        atomic_number = atomic_numbers[symbol]
    else:
        atomic_number = symbol

    # Interpret noshells
    if noshells < 1:
        raise ValueError("The number of shells must be equal to or greater than one.")

    # Interpret lattice constant
    if latticeconstant is None:
        if reference_states[atomic_number]['symmetry'] in ['fcc', 'bcc', 'sc']:
            lattice_constant = reference_states[atomic_number]['a']
        else:
            raise NotImplementedError(("Cannot guess lattice constant of a %s element." %
                                      (reference_states[atomic_number]['symmetry'],)))
    else:
        if isinstance(latticeconstant, (int, float)):
            lattice_constant = latticeconstant
        else:
            raise ValueError("Lattice constant must be of type int or float.")

    t = 0.5 + np.sqrt(5)/2.0

    verticies = np.array([[t, 0., 1.],
                          [t, 0., -1.],
                          [-t, 0., 1.],
                          [-t, 0., -1.],
                          [1., t, 0.],
                          [-1., t, 0.],
                          [1., -t, 0.],
                          [-1., -t, 0.],
                          [0., 1., t],
                          [0., -1., t],
                          [0., 1., -t],
                          [0., -1., -t]])

    positions = []
    tags = []
    positions.append(np.zeros(3))
    tags.append(1)

    for n in range(1, noshells):
        #Construct square edges (6)
        for k in range(0, 12, 2):
            v1 = verticies[k]
            v2 = verticies[k+1]
            for i in range(n+1):
                pos = i*v1 + (n-i)*v2
                positions.append(pos)
                tags.append(n + 1)

        #Construct triangle planes (12)
        if n > 1:
            map = {0: (8, 9), 1: (10, 11),
                   2: (8, 9), 3: (10, 11),
                   4: (0, 1), 5: (2, 3),
                   6: (0, 1), 7: (2, 3),
                   8: (4, 5), 9: (6, 7),
                   10: (4, 5), 11: (6, 7)}

            for k in range(0, 12):
                v0 = n*verticies[k]
                v1 = (verticies[map[k][0]] - verticies[k])
                v2 = (verticies[map[k][1]] - verticies[k])
                for i in range(n):
                    for j in range(n-i):
                        if i == 0 and j == 0:
                            continue
                        pos = v0 + i*v1 + j*v2
                        positions.append(pos)
                        tags.append(n + 1)

        #Fill missing triangle planes (8)
        if n > 2:
            map = {0: (9, 6, 8, 4,),
                   1: (11, 6, 10, 4),
                   2: (9, 7, 8, 5,),
                   3: (11, 7, 10, 5)}

            for k in range(0, 4):
                v0 = n*verticies[k]
                v1 = (verticies[map[k][0]] - verticies[k])
                v2 = (verticies[map[k][1]] - verticies[k])
                v3 = (verticies[map[k][2]] - verticies[k])
                v4 = (verticies[map[k][3]] - verticies[k])
                for i in range(1, n):
                    for j in range(1, n-i):
                        pos = v0 + i*v1 + j*v2
                        positions.append(pos)
                        tags.append(n + 1)
                        pos = v0 + i*v3 + j*v4
                        positions.append(pos)
                        tags.append(n + 1)

    # Scale the positions
    scaling_factor = lattice_constant / np.sqrt(2*(1 + t**2))
    positions = np.array(positions) * scaling_factor

    # Fit the cell, so it only just consist the atoms
    min = positions.min(axis=0)
    max = positions.max(axis=0)
    cell = max - min
    positions = positions - min

    symbols = [atomic_number] * len(positions)
    return Atoms(symbols=symbols, positions=positions, tags=tags, cell=cell)

