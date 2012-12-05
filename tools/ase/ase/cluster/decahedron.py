import numpy as np

from ase import Atoms
from ase.data import atomic_numbers, reference_states

def Decahedron(symbol, p, q, r, latticeconstant=None):
    """
    Returns a cluster in the decahedra class.

    Prameters
    ---------
    symbol: Chemical symbol (or atomic number) of the element.

    p: Number of atoms on the (100) facets perpendicular to the five
    fold axis.
    
    q: Number of atoms on the (100) facets parallel to the five fold
    axis. q = 1 corresponds to no visible (100) facets.

    r: Depth of the Marks re-entrence at the pentagon corners. r = 0
    corresponds to no re-entrence.

    latticeconstant (optional): The lattice constant. If not given,
    then it is extracted form ase.data.
    """

    # Interpret symbol
    if isinstance(symbol, str):
        atomic_number = atomic_numbers[symbol]
    else:
        atomic_number = symbol

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

    # Check values of p, q, r
    if p < 1 or q < 1:
        raise ValueError("p and q must be greater than 0.")

    if r < 0:
        raise ValueError("r must be greater than or equal to 0.")

    # Defining constants
    t = 2.0*np.pi/5.0
    b = lattice_constant/np.sqrt(2.0)
    a = b*np.sqrt(3.0)/2.0

    verticies = a * np.array([[np.cos(np.pi/2.), np.sin(np.pi/2.), 0.],
                              [np.cos(t*1. + np.pi/2.), np.sin(t*1. + np.pi/2.), 0.],
                              [np.cos(t*2. + np.pi/2.), np.sin(t*2. + np.pi/2.), 0.],
                              [np.cos(t*3. + np.pi/2.), np.sin(t*3. + np.pi/2.), 0.],
                              [np.cos(t*4. + np.pi/2.), np.sin(t*4. + np.pi/2.), 0.]])

    # Number of atoms on the five fold axis and a nice constant
    h = p + q + 2*r - 1
    g = h - q + 1 # p + 2*r

    positions = []
    # Make the five fold axis
    for j in range(h):
        pos = np.array([0.0, 0.0, j*b - (h-1)*b/2.0])
        positions.append(pos)

    # Make pentagon rings around the five fold axis
    for n in range(1, h):
        # Condition for (100)-planes
        if n < g:
            for m in range(5):
                v1 = verticies[m-1]
                v2 = verticies[m]
                for i in range(n):
                    # Condition for marks re-entrence
                    if n - i < g - r and i < g - r:
                        for j in range(h-n):
                            pos = (n-i)*v1 + i*v2
                            pos += np.array([0.0, 0.0, j*b - (h-n-1)*b/2.0])
                            positions.append(pos)

    # Fit the cell, so it only just consist the atoms
    min = np.zeros(3)
    max = np.zeros(3)
    axes = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    for i in range(3):
        r = np.dot(positions, axes[i])
        min[i] = r.min()
        max[i] = r.max()
    cell = max - min
    positions = np.array(positions) - min

    symbols = [atomic_number] * len(positions)
    return Atoms(symbols=symbols, positions=positions, cell=cell)
