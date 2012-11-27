from math import pi, cos, sin, sqrt, acos

from ase.atoms import Atoms
from ase.parallel import paropen

def read_sdf(fileobj):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    lines = fileobj.readlines()
    # first three lines header
    del lines[:3]

    # 
    L1 = lines.pop(0).split()
    natoms = int(L1[0])
    positions = []
    symbols = []
    for line in lines[:natoms]:
        x, y, z, symbol = line.split()[:4]
        symbols.append(symbol)
        positions.append([float(x), float(y), float(z)])
    return Atoms(symbols=symbols, positions=positions)
