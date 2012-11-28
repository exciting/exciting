from math import pi, cos, sin, sqrt, acos
import numpy as np

from ase.data import chemical_symbols
from ase.atoms import Atoms
from ase.parallel import paropen

iwm_symbols = {'1' : 'C',
               '2' : 'Au',
               '5' : 'Ag'}

def read_iwm(fileobj, index=-1):
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    lines = fileobj.readlines()
    L1 = lines[1].split()
    if len(L1) == 1:
        del lines[:3]
        natoms = int(L1[0])
    else:
        natoms = len(lines)
    images = []

    positions = []
    symbols = []
    for line in lines[:natoms]:
        symbol, mass, x, y, z = line.split()[:5]
        if symbol in iwm_symbols:
            symbols.append(iwm_symbols[symbol])
        else:
            symbols.append(chemical_symbols[int(symbol)])
        positions.append([float(x), float(y), float(z)])
    
    del(lines[natoms:3 * natoms + 3])
    
    cell = []
    for line in lines[natoms:natoms+3]:
        x, y, z = line.split()[:3]
        cell.append(np.array([float(x), float(y), float(z)]))
      
    images.append(Atoms(symbols=symbols, positions=positions, cell=cell))

    return images[index]

