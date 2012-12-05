"""
This module contains functionality for reading an ASE
Atoms object in V_Sim ascii format.

"""

import os

def read_v_sim(filename='demo.ascii'):
    """Import V_Sim input file.

    Reads cell, atom positions, etc. from v_sim ascii file
    """

    from ase import Atoms, units
    from ase.lattice.spacegroup import cell
    import re

    if isinstance(filename, str):
        f = open(filename)
    else: # Assume it's a file-like object
        f = filename

    comment = f.readline()

    line = f.readline() + ' ' + f.readline()
    box = line.split()
    for i in range(len(box)):
        box[i] = float(box[i])

    keywords = []
    positions= []
    symbols  = []
    unit     = 1.0

    re_comment = re.compile('^\s*[#!]')
    re_node    = re.compile('^\s*\S+\s+\S+\s+\S+\s+\S+')

    while(True):
        line = f.readline()

        if line == '':
            break # EOF

        p = re_comment.match(line)
        if p != None:
            # remove comment character at the beginning of line
            line = line[p.end():].replace(',', ' ').lower()
            if line[:8] == "keyword:":
                keywords.extend(line[8:].split())

        elif(re_node.match(line)):
            unit = 1.0
            if not ("reduced" in keywords) :
                if ("bohr" in keywords) or ("bohrd0" in keywords) or ("atomic" in keywords) or ("atomicd0" in keywords):
                    unit = units.Bohr

            fields = line.split()
            positions.append([unit*float(fields[0]),
                              unit*float(fields[1]),
                              unit*float(fields[2])])
            symbols.append(fields[3])

    f.close()

    if ("surface" in keywords) or ("freeBC" in keywords):
        raise NotImplementedError

    # create atoms object based on the information
    if ("angdeg" in keywords) :
        cell = cellpar_to_cell(box)
    else:
        unit = 1.0
        if ("bohr" in keywords) or ("bohrd0" in keywords) or ("atomic" in keywords) or ("atomicd0" in keywords):
            unit = units.Bohr
        cell = [[unit*box[0],         0.0,         0.0],
                [unit*box[1], unit*box[2],         0.0],
                [unit*box[3], unit*box[4], unit*box[5]]]

    if ("reduced" in keywords) :
        atoms = Atoms(cell=cell, scaled_positions=positions)
    else :
        atoms = Atoms(cell=cell, positions=positions)

    atoms.set_chemical_symbols(symbols)

    return atoms
