"""
IO support for the Gaussian cube format.

See the format specifications on:
http://local.wasp.uwa.edu.au/~pbourke/dataformats/cube/
"""


import numpy as np

from ase.atoms import Atoms
from ase.units import Bohr
from ase.parallel import paropen


def write_cube(fileobj, atoms, data=None):
    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')
        
    if isinstance(atoms, list):
        if len(atoms) > 1:
            raise ValueError('Can only write one configuration '
                             'to a cube file!')
        atoms = atoms[0]

    if data is None:
        data = np.ones((2, 2, 2))
    data = np.asarray(data)

    if data.dtype == complex:
        data = np.abs(data)

    fileobj.write('cube file from ase\n')
    fileobj.write('OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n')

    cell = atoms.get_cell()
    shape = np.array(data.shape)

    corner = np.zeros(3)
    for i in range(3):
        if shape[i] % 2 == 1:
            shape[i] += 1
            corner += cell[i] / shape[i] / Bohr

    fileobj.write('%5d%12.6f%12.6f%12.6f\n' % (len(atoms), corner[0],
                                               corner[1], corner[2]))

    for i in range(3):
        n = data.shape[i]
        d = cell[i] / shape[i] / Bohr
        fileobj.write('%5d%12.6f%12.6f%12.6f\n' % (n, d[0], d[1], d[2]))

    positions = atoms.get_positions() / Bohr
    numbers = atoms.get_atomic_numbers()
    for Z, (x, y, z) in zip(numbers, positions):
        fileobj.write('%5d%12.6f%12.6f%12.6f%12.6f\n' % (Z, 0.0, x, y, z)) 

    data.tofile(fileobj, sep='\n', format='%e')


def read_cube(fileobj, index=-1, read_data=False):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    readline = fileobj.readline
    readline()
    axes = ['XYZ'.index(s[0]) for s in readline().split()[2::3]]
    if axes == []:
        axes = [0, 1, 2]
    line = readline().split()
    natoms = int(line[0])
    corner = [Bohr * float(x) for x in line[1:]]

    cell = np.empty((3, 3))
    shape = []
    for i in range(3):
        n, x, y, z = [float(s) for s in readline().split()]
        shape.append(n)
        if n % 2 == 1:
            n += 1
        cell[i] = n * Bohr * np.array([x, y, z])
        
    numbers = np.empty(natoms, int)
    positions = np.empty((natoms, 3))
    for i in range(natoms):
        line = readline().split()
        numbers[i] = int(line[0])
        positions[i] = [float(s) for s in line[2:]]

    positions *= Bohr
    atoms = Atoms(numbers=numbers, positions=positions, cell=cell)

    if read_data:
        data = np.array([float(s)
                         for s in fileobj.read().split()]).reshape(shape)
        if axes != [0, 1, 2]:
            data = data.transpose(axes).copy()
        return data, atoms

    return atoms


def read_cube_data(fileobj):
    return read_cube(fileobj, read_data=True)
