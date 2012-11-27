import numpy as np 

from ase.atoms import Atom, Atoms
from ase.parallel import paropen

"""Module to read and write atoms in PDB file format"""


def read_pdb(fileobj, index=-1):
    """Read PDB files.

    The format is assumed to follow the description given in
    http://www.wwpdb.org/documentation/format32/sect9.html."""
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    images = []
    atoms = Atoms()
    for line in fileobj.readlines():
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                # Atom name is arbitrary and does not necessarily contain the element symbol.
                # The specification requires the element symbol to be in columns 77+78.
                symbol = line[76:78].strip().lower().capitalize()
                words = line[30:55].split()
                position = np.array([float(words[0]), 
                                     float(words[1]),
                                     float(words[2])])
                atoms.append(Atom(symbol, position))
            except:
                pass
        if line.startswith('ENDMDL'):
            images.append(atoms)
            atoms = Atoms()
    if len(images) == 0:
        images.append(atoms)
    return images[index]

def write_pdb(fileobj, images):
    """Write images to PDB-file.

    The format is assumed to follow the description given in
    http://www.wwpdb.org/documentation/format32/sect9.html."""
    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    if images[0].get_pbc().any():
        from ase.lattice.spacegroup.cell import cell_to_cellpar
        cellpar = cell_to_cellpar( images[0].get_cell())
        # ignoring Z-value, using P1 since we have all atoms defined explicitly
        format = 'CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1\n'
        fileobj.write(format % (cellpar[0], cellpar[1], cellpar[2], cellpar[3], cellpar[4], cellpar[5]))

    #         1234567 123 6789012345678901   89   67   456789012345678901234567 890
    format = 'ATOM  %5d %4s MOL     1    %8.3f%8.3f%8.3f  1.00  0.00          %2s  \n'

    # RasMol complains if the atom index exceeds 100000. There might
    # be a limit of 5 digit numbers in this field.
    MAXNUM = 100000

    symbols = images[0].get_chemical_symbols()
    natoms = len(symbols)
    
    for n,atoms in enumerate(images):
        fileobj.write('MODEL     '+str(n+1)+'\n')
        p = atoms.get_positions()
        for a in range(natoms):
            x, y, z = p[a]
            fileobj.write(format % (a % MAXNUM, symbols[a], x, y, z, symbols[a].rjust(2)))
        fileobj.write('ENDMDL\n')
