"""Extension to ASE: read and write structures in GEN format

Refer to DFTB+ manual for GEN format description.

Note: GEN format only supports single snapshot.
"""

from ase.atoms import Atoms
from ase.parallel import paropen


def read_gen(fileobj):
    """Read structure in GEN format (refer to DFTB+ manual).
       Multiple snapshot are not allowed. """
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    image = Atoms()    
    lines = fileobj.readlines()
    line = lines[0].split()
    natoms = int(line[0])
    if line[1] == 'S':
        supercell = True
    elif line[1] == 'C':
        supercell = False
    else:
        raise IOError('Error in line #1: only C (Cluster) or S (Supercell) ' +
                      'are valid options')

    # Read atomic symbols
    line = lines[1].split()
    # Define a dictionary with symbols-id
    symboldict = dict()
    symbolid = 1
    for symb in line:
        symboldict[symbolid] = symb
        symbolid += 1

    # Read atoms (GEN format supports only single snapshot)
    del lines[:2]
    positions = []
    symbols = []
    for line in lines[:natoms]:
        dummy, symbolid, x, y, z = line.split()[:5]
        symbols.append(symboldict[int(symbolid)])
        positions.append([float(x), float(y), float(z)])
    image = Atoms(symbols=symbols, positions=positions)
    del lines[:natoms]

    # If Supercell, parse periodic vectors
    if not supercell:
        return image
    else:
        # Dummy line: line after atom positions is not uniquely defined 
        # in gen implementations, and not necessary in DFTB package
        del lines[:1]
        image.set_pbc([True, True, True])
        p = []
        for i in range(3):
            x, y, z = lines[i].split()[:3]
            p.append([float(x), float(y), float(z)])
        image.set_cell([(p[0][0], p[0][1], p[0][2]), (p[1][0], p[1][1], 
                p[1][2]), (p[2][0], p[2][1], p[2][2])])
        return image
        

def write_gen(fileobj, images):
    """Write structure in GEN format (refer to DFTB+ manual).
       Multiple snapshots are not allowed. """
    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    # Images is kept in a list but a size > 0 is not allowed
    # as GEN format doesn't support multiple snapshots.
    # Images is used as a list for consistency with the other 
    # output modules
    if len(images) != 1:
        raise ValueError('images contains more than one structure\n' +
                         'GEN format supports only single snapshot output')

    symbols = images[0].get_chemical_symbols()

    # Define a dictionary with symbols-id
    symboldict = dict() 
    for sym in symbols:
        if not (sym in symboldict):
            symboldict[sym] = len(symboldict) + 1
    # An ordered symbol list is needed as ordered dictionary
    # is just available in python 2.7
    orderedsymbols = list(['null'] * len(symboldict.keys()))
    for sym in symboldict.keys():
        orderedsymbols[symboldict[sym] - 1] = sym
    
    
    # Check whether the structure is periodic
    # GEN cannot describe periodicity in one or two direction,
    # a periodic structure is considered periodic in all the
    # directions. If your structure is not periodical in all
    # the directions, be sure you have set big periodicity
    # vectors in the non-periodic directions
    if images[0].pbc.any():
        pb_flag = 'S'
    else:
        pb_flag = 'C'

    natoms = len(symbols)
    ind = 0
    for atoms in images:
        fileobj.write('%d  %-5s\n' % (natoms, pb_flag))
        for s in orderedsymbols:
            fileobj.write('%-5s' % s)
        fileobj.write('\n')
        for sym, (x, y, z) in zip(symbols, atoms.get_positions()):
            ind += 1
            symbolid = symboldict[sym]
            fileobj.write('%-6d %d %22.15f %22.15f %22.15f\n' % (ind, 
                          symbolid, x, y, z))
    if images[0].pbc.any():
        fileobj.write('%22.15f %22.15f %22.15f \n' % (0.0, 0.0, 0.0))
        fileobj.write('%22.15f %22.15f %22.15f \n' % 
                      (images[0].get_cell()[0][0], 
                       images[0].get_cell()[0][1], 
                       images[0].get_cell()[0][2]))
        fileobj.write('%22.15f %22.15f %22.15f \n' % 
                      (images[0].get_cell()[1][0], 
                       images[0].get_cell()[1][1], 
                       images[0].get_cell()[1][2]))
        fileobj.write('%22.15f %22.15f %22.15f \n' % 
                      (images[0].get_cell()[2][0], 
                       images[0].get_cell()[2][1], 
                       images[0].get_cell()[2][2]))
