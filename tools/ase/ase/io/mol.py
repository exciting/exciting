from math import pi, cos, sin, sqrt, acos

from ase.atoms import Atoms
from ase.parallel import paropen


def read_mol(fileobj, index=-1):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    lines = fileobj.readlines()
    del(lines[:3])
    L1 = lines[0].split()
    del(lines[0])
    natoms = int(L1[0])
    positions = []
    symbols = []
    for line in lines[:natoms]:
            x, y, z, symbol = line.split()[:4]
            symbols.append(symbol)
            positions.append([float(x), float(y), float(z)])
    return Atoms(symbols=symbols, positions=positions)

