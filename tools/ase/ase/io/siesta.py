from numpy import zeros
from os import fstat
from re import compile

from ase.io.fortranfile import FortranFile


def read_rho(fname):
    "Read unformatted Siesta charge density file"

    # TODO:
    # 
    # Handle formatted and NetCDF files.
    #
    # Siesta source code (at least 2.0.2) can possibly also 
    # save RHO as a _formatted_ file (the source code seems
    # prepared, but there seems to be no fdf-options for it though).
    # Siesta >= 3 has support for saving RHO as a NetCDF file
    # (according to manual)

    fh = FortranFile(fname)
    
    # Read (but ignore) unit cell vectors
    x = fh.readReals('d')
    if len(x) != 3 * 3: 
        raise IOError('Failed to read cell vectors')
        
    # Read number of grid points and spin components
    x = fh.readInts()
    if len(x) != 4:
        raise IOError('Failed to read grid size')
    gpts = x  # number of 'X', 'Y', 'Z', 'spin' gridpoints 
    
    rho = zeros(gpts)
    for ispin in range(gpts[3]):
        for n3 in range(gpts[2]):
            for n2 in range(gpts[1]):
                x = fh.readReals('f')
                if len(x) != gpts[0]:
                    raise IOError('Failed to read RHO[:,%i,%i,%i]' %
                                  (n2, n3, ispin))
                rho[:, n2, n3, ispin] = x
    
    fh.close()
    return rho



#
# Helper functions for read_fdf
#
_label_strip_re = compile(r'[\s._-]')
def _labelize(raw_label):
    # Labels are case insensitive and -_. should be ignored, lower and strip it
    return _label_strip_re.sub('', raw_label).lower()

def _is_block(val):
    # Tell whether value is a block-value or an ordinary value.
    # A block is represented as a list of lists of strings,
    # and a ordinary value is represented as a list of strings
    if type(val) is list and \
       len(val) > 0 and \
       type(val[0]) is list:
        return True
    return False
    
def _get_stripped_lines(fd):
    # Remove comments, leading blanks, and empty lines
    return filter(None, [L.split('#')[0].strip() for L in fd])

def _read_fdf_lines(file, inodes=[]):
    # Read lines and resolve includes

    if type(file) is str:
        file = open(file, 'r')
    fst = fstat(file.fileno())
    inode = (fst.st_dev, fst.st_ino)
    if inode in inodes:
        raise IOError('Cyclic include in fdf file')
    inodes = inodes + [inode]

    lbz = _labelize
    
    lines = []
    for L in _get_stripped_lines(file):
        w0 = lbz(L.split(None, 1)[0])

        if w0 == '%include':
            # Include the contents of fname
            fname = L.split(None, 1)[1].strip()
            lines += _read_fdf_lines(fname, inodes)

        elif '<' in L:
            L, fname = L.split('<', 1)
            w = L.split()
            fname = fname.strip()

            if w0 == '%block':
                # "%block label < filename" means that the block contents should be read from filename 
                if len(w) != 2:
                    raise IOError('Bad %%block-statement "%s < %s"' % (L, fname))
                label = lbz(w[1])
                lines.append('%%block %s' % label)
                lines += _get_stripped_lines(open(fname))
                lines.append('%%endblock %s' % label)
            else:
                # "label < filename.fdf" means that the label (_only_ that label) is to be resolved from filename.fdf
                label = lbz(w[0])
                fdf = _read_fdf(fname, inodes)
                if label in fdf:
                    if _is_block(fdf[label]):
                        lines.append('%%block %s' % label)
                        lines += [' '.join(x) for x in fdf[label]]
                        lines.append('%%endblock %s' % label)
                    else:
                        lines.append('%s %s' % (label, ' '.join(fdf[label])))
                #else: label unresolved! One should possibly issue a warning about this!
        else:
            # Simple include line L
            lines.append(L)
    return lines

#
# The reason for creating a separate _read_fdf is simply to hide the inodes-argument
#
def _read_fdf(fname, inodes=[]):
    # inodes is used to detect cyclic includes
    fdf = {}
    lbz = _labelize
    lines = _read_fdf_lines(fname, inodes)
    while lines:
        w = lines.pop(0).split(None, 1)
        if lbz(w[0]) == '%block':
            # Block value
            if len(w) == 2:
                label = lbz(w[1])
                content = []
                while True:
                    if len(lines) == 0:
                        raise IOError('Unexpected EOF reached in %s, '
                                      'un-ended block %s' % (fname, label))                            
                    w = lines.pop(0).split()
                    if lbz(w[0]) == '%endblock' and lbz(w[1]) == label:
                        break
                    content.append(w)
                        
                if not label in fdf:
                    # Only first appearance of label is to be used
                    fdf[label] = content
            else:
                raise IOError('%%block statement without label' )
        else:
            # Ordinary value
            label = lbz(w[0])
            if len(w) == 1:
                # Siesta interpret blanks as True for logical variables
                fdf[label] = []
            else:
                fdf[label] = w[1].split()
    return fdf


def read_fdf(fname):
    """Read a siesta style fdf-file.

    The data is returned as a dictionary
    ( label:value ).
    
    All labels are converted to lower case characters and
    are stripped of any '-', '_', or '.'.
    
    Ordinary values are stored as a list of strings (splitted on WS),
    and block values are stored as list of lists of strings
    (splitted per line, and on WS).
    If a label occurres more than once, the first occurrence
    takes precedence.

    The implementation applies no intelligence, and does not
    "understand" the data or the concept of units etc.
    Values are never parsed in any way, just stored as
    split strings.
    
    The implementation tries to comply with the fdf-format
    specification as presented in the siesta 2.0.2 manual.

    An fdf-dictionary could e.g. look like this::

        {'atomiccoordinatesandatomicspecies': [
              ['4.9999998', '5.7632392', '5.6095972', '1'],
              ['5.0000000', '6.5518100', '4.9929091', '2'],
              ['5.0000000', '4.9746683', '4.9929095', '2']],
         'atomiccoordinatesformat': ['Ang'],
         'chemicalspecieslabel': [['1', '8', 'O'],
                                  ['2', '1', 'H']],
         'dmmixingweight': ['0.1'],
         'dmnumberpulay': ['5'],
         'dmusesavedm': ['True'],
         'latticeconstant': ['1.000000', 'Ang'],
         'latticevectors': [
              ['10.00000000', '0.00000000', '0.00000000'],
              ['0.00000000', '11.52647800', '0.00000000'],
              ['0.00000000', '0.00000000', '10.59630900']],
         'maxscfiterations': ['120'],
         'meshcutoff': ['2721.139566', 'eV'],
         'numberofatoms': ['3'],
         'numberofspecies': ['2'],
         'paobasissize': ['dz'],
         'solutionmethod': ['diagon'],
         'systemlabel': ['H2O'],
         'wavefunckpoints': [['0.0', '0.0', '0.0']],
         'writedenchar': ['T'],
         'xcauthors': ['PBE'],
         'xcfunctional': ['GGA']}

    """

    return _read_fdf(fname)


def read_struct(fname):
    """Read a siesta struct file"""
    from ase.atoms import Atoms, Atom

    f = open(fname, 'r')

    cell = []
    for i in range(3):
        cell.append([float(x) for x in f.readline().split()])

    natoms = int(f.readline())

    atoms = Atoms()
    for atom in f:
        Z, pos_x, pos_y, pos_z = atom.split()[1:]
        atoms.append(Atom(int(Z), position = (float(pos_x), float(pos_y), float(pos_z))))

    if len(atoms) != natoms:
        raise IOError('Badly structured input file')
    
    atoms.set_cell(cell, scale_atoms = True)

    return atoms
    
