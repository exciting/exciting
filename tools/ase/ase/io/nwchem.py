from cStringIO import StringIO
from ase.atoms import Atoms
from ase.io.xyz import read_xyz

def read_nwchem(filename):
    """Method to read geometry from a nwchem output
    """
    from ase import Atoms, Atom

    f = filename
    if isinstance(filename, str):
        f = open(filename)

    lines = f.readlines()

    i = 0
    while i < len(lines):
        if lines[i].find('XYZ format geometry') >=0:
            natoms = int(lines[i + 2].split()[0])
            string = ''
            for j in range(2, natoms + 4):
                xyzstring = lines[i + j]
                symbol = xyzstring.split()[0].strip()
                # replace bq ghost with X: MDTMP can we do better?
                if symbol.startswith('bq'):
                    xyzstring = xyzstring.replace(symbol, 'X')
                string += xyzstring
            atoms = read_xyz(StringIO(string))
            i += natoms + 4
        else:
            i += 1

    if type(filename) == str:
        f.close()

    return atoms

def write_nwchem(filename, atoms, geometry=None):
    """Method to write nwchem coord file
    """

    import numpy as np

    if isinstance(filename, str):
        f = open(filename, 'w')
    else: # Assume it's a 'file-like object'
        f = filename

    # autosym and autoz are defaults
    # http://www.nwchem-sw.org/index.php/Geometry
    # geometry noautoz results in higher memory demand!
    # http://www.emsl.pnl.gov/docs/nwchem/nwchem-support/2010/10/0060.RE:_NWCHEM_Geometry_problem_fwd_
    if geometry is not None:
        f.write('geometry ' + str(geometry) + '\n')
    else:
        f.write('geometry\n')
    for atom in atoms:
        if atom.tag == -71: # 71 is ascii G (Ghost)
            symbol = 'bq' + atom.symbol
        else:
            symbol = atom.symbol
        f.write('  ' + symbol + ' ' +
                str(atom.position[0]) + ' ' +
                str(atom.position[1]) + ' ' +
                str(atom.position[2]) + '\n' )
    f.write('end\n')
