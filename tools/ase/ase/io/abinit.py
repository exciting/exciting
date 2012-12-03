"""
This module contains functionality for reading an ASE
Atoms object in ABINIT input format.

"""

import os

def read_abinit(filename='abinit.in'):
    """Import ABINIT input file.

    Reads cell, atom positions, etc. from abinit input file
    """

    from ase import Atoms, units

    if isinstance(filename, str):
        f = open(filename)
    else: # Assume it's a file-like object
        f = filename

    lines = f.readlines()
    if type(filename) == str:
        f.close()

    full_file = ''
    for line in lines:
        if '#' in line:
            meat, comment = line.split('#')
        else:
            meat = line
        full_file = full_file + meat + ' '

    full_file.strip()
    tokens = full_file.lower().split()

    # note that the file can not be scanned sequentially

    index = tokens.index("acell")
    unit = 1.0
    if(tokens[index+4].lower()[:3] != 'ang'):
        unit = units.Bohr
    acell = [unit*float(tokens[index+1]),
             unit*float(tokens[index+2]),
             unit*float(tokens[index+3])]

    index = tokens.index("natom")
    natom = int(tokens[index+1])

    index = tokens.index("ntypat")
    ntypat = int(tokens[index+1])

    index = tokens.index("typat")
    typat = []
    for i in range(natom):
        typat.append(int(tokens[index+1+i]))

    index = tokens.index("znucl")
    znucl = []
    for i in range(ntypat):
        znucl.append(int(tokens[index+1+i]))

    index = tokens.index("rprim")
    rprim = []
    for i in range(3):
        rprim.append([acell[i]*float(tokens[index+3*i+1]),
                      acell[i]*float(tokens[index+3*i+2]),
                      acell[i]*float(tokens[index+3*i+3])])

    # create a list with the atomic numbers
    numbers = []
    for i in range(natom):
        ii = typat[i] - 1
        numbers.append(znucl[ii])

    # now the positions of the atoms
    if "xred" in tokens:
        index = tokens.index("xred")
        xred = []
        for i in range(natom):
            xred.append([float(tokens[index+3*i+1]),
                         float(tokens[index+3*i+2]),
                         float(tokens[index+3*i+3])])
        atoms = Atoms(cell=rprim, scaled_positions=xred, numbers=numbers)
        return atoms

    index = None
    if "xcart" in tokens:
        index = tokens.index("xcart")
        unit = units.Bohr
    elif "xangs" in tokens:
        unit = 1.0
        index = tokens.index("xangs")

    if(index != None):
        xangs = []
        for i in range(natom):
            xangs.append([unit*float(tokens[index+3*i+1]),
                          unit*float(tokens[index+3*i+2]),
                          unit*float(tokens[index+3*i+3])])
        atoms = Atoms(cell=rprim, positions=xangs, numbers=numbers)
        return atoms

    raise IOError("No xred, xcart, or xangs keyword in abinit input file")


def write_abinit(filename, atoms, cartesian=False, long_format=True):
    """Method to write abinit input files."""

    import numpy as np
    from ase import Atoms, data

    if isinstance(filename, str):
        f = open(filename, 'w')
    else: # Assume it's a 'file-like object'
        f = filename

    if isinstance(atoms, (list, tuple)):
        if len(atoms) > 1:
            raise RuntimeError("Don't know how to save more than "+
                            "one image to input")
        else:
            atoms = atoms[0]

    # Write atom positions in scaled or cartesian coordinates
    if cartesian:
        coord = atoms.get_positions()
    else:
        coord = atoms.get_scaled_positions()

    # let us order the atoms according to chemical symbol
    ind = np.argsort(atoms.get_chemical_symbols())
    symbols = np.array(atoms.get_chemical_symbols())[ind]
    coord = coord[ind]

    # and now we count how many atoms of which type we have
    sc = []
    psym = symbols[0]
    count = 0
    for sym in symbols:
        if sym != psym:
            sc.append((psym, count))
            psym = sym
            count = 1
        else:
            count += 1
    sc.append((psym, count))

    f.write('\n# Definition of the atom types\n')
    f.write("ntypat  " + str(len(sc)) + "\n")
    f.write("znucl  ")
    for specie in sc:
        f.write(str(data.atomic_numbers[specie[0]]) + " ")
    f.write('\n')

    f.write('\n# Definition of the atoms\n')
    f.write('natom  ' + str(len(symbols)) + '\n')
    f.write('typat  ')
    typat = 1
    for specie in sc:
        for natom in range(specie[1]):
            f.write(str(typat) + ' ')
        typat = typat + 1
    f.write('\n')

    f.write('\n# Definition of the unit cell\n')
    f.write('acell\n')
    f.write('%.14f %.14f %.14f Angstrom\n' %  (1.0, 1.0, 1.0))
    f.write('\n')
    f.write('rprim\n')
    if long_format:
        latt_form = ' %21.16f'
    else:
        latt_form = ' %11.6f'

    for vec in atoms.get_cell():
        f.write(' ')
        for el in vec:
            f.write(latt_form % el)
        f.write('\n')
    f.write('\n')

    # Write atom positions in scaled or cartesian coordinates
    if cartesian:
        f.write('xangst\n')
    else:
        f.write('xred\n')

    if long_format:
        cform = ' %19.16f'
    else:
        cform = ' %9.6f'

    for iatom, atom in enumerate(coord):
        f.write(' ')
        for dcoord in atom:
            f.write(cform % dcoord)
        f.write('\n')

    if type(filename) == str:
        f.close()
