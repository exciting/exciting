# Copyright (C) 2012, Jesper Friis
# (see accompanying license files for ASE).
"""
Determines space group of an atoms object using the FINDSYM program
from the ISOTROPY (http://stokes.byu.edu/iso/isotropy.html) software
package by H. T. Stokes and D. M. Hatch, Brigham Young University,
USA.

In order to use this module, you have to download the ISOTROPY package
from http://stokes.byu.edu/iso/isotropy.html and set the environment
variable ISODATA to the path of the directory containing findsym
and data_space.txt (NB: the path should end with a slash (/)).


Example
-------
>>> from ase.lattice.spacegroup import crystal
>>> from ase.utils.geometry import cut

# Start with simple fcc Al
>>> al = crystal('Al', [(0,0,0)], spacegroup=225, cellpar=4.05)
>>> d = findsym(al)
>>> d['spacegroup']
225

# No problem with a more complex structure...
>>> skutterudite = crystal(('Co', 'Sb'),
...                        basis=[(0.25,0.25,0.25), (0.0, 0.335, 0.158)],
...                        spacegroup=204, 
...                        cellpar=9.04)
>>> d = findsym(skutterudite)
>>> d['spacegroup']
204

# ... or a non-conventional cut
slab = cut(skutterudite, a=(1, 1, 0), b=(0, 2, 0), c=(0, 0, 1))
d = findsym(slab)
>>> d['spacegroup']
204
"""

import os
import subprocess

import numpy as np
import ase

__all__ = ['findsym', 'unique']


def make_input(atoms, tol=1e-3, centering='P', types=None):
    """Returns input to findsym.  See findsym() for a description of
    the arguments."""
    if types is None:
        types = atoms.numbers
    s = []
    s.append(atoms.get_name())
    s.append('%g  tolerance' % tol)
    s.append('2    form of lattice parameters: to be entered as lengths '
             'and angles')
    s.append('%g %g %g %g %g %g    a,b,c,alpha,beta,gamma' % 
             tuple(ase.lattice.spacegroup.cell.cell_to_cellpar(atoms.cell)))
    s.append('2    form of vectors defining unit cell')  # ??
    s.append('%s    centering (P=unknown)' % centering)
    s.append('%d   number of atoms in primitive unit cell' % len(atoms))
    s.append(' '.join(str(n) for n in types) + '   type of each atom')
    for p in atoms.get_scaled_positions():
        s.append('%10.5f  %10.5f  %10.5f' % tuple(p))
    return '\n'.join(s)


def run(atoms, tol=1e-3, centering='P', types=None, isodata_dir=None):
    """Runs FINDSYM and returns its standard output."""
    if isodata_dir is None:
        isodata_dir = os.getenv('ISODATA')
    if isodata_dir is None:
        isodata_dir = '.'
    isodata_dir = os.path.normpath(isodata_dir)
    findsym = os.path.join(isodata_dir, 'findsym')
    data_space = os.path.join(isodata_dir, 'data_space.txt')
    for path in findsym, data_space:
        if not os.path.exists(path):
            raise IOError('no such file: %s. Have you set the ISODATA '
                          'environment variable to the directory containing '
                          'findsym and data_space.txt?' % path)
    env = os.environ.copy()
    env['ISODATA'] = isodata_dir + os.sep
    p = subprocess.Popen([findsym], stdin=subprocess.PIPE, 
                         stdout=subprocess.PIPE, env=env)
    stdout = p.communicate(make_input(atoms, tol, centering, types))[0]
    #if os.path.exists('findsym.log'):
    #    os.remove('findsym.log')
    return stdout
    

def parse(output):
    """Parse output from FINDSYM (Version 3.2.3, August 2007) and
    return a dict.  See docstring for findsym() for a description of
    the tokens."""
    import StringIO
    d = {}
    f = StringIO.StringIO(output)

    line = f.readline()
    while not line.startswith('Lattice parameters'):
        line = f.readline()
        assert line, 'Cannot find line starting with "Lattice parameters"'
    d['cellpar'] = np.array([float(v) for v in f.readline().split()])
    #d['centering'] = f.readline().split()[1]  # seems not to make sence...

    # Determine number of atoms from atom types, since the number of
    # atoms is written with only 3 digits, which crashes the parser
    # for more than 999 atoms
    while not line.startswith('Type of each atom'):
        line = f.readline()
        assert line, 'Cannot find line starting with "Type of each atom"'
    line = f.readline()
    types = []
    while not line.startswith('Position of each atom'):
        types.extend([int(n) for n in line.split()])
        line = f.readline()
        assert line, 'Cannot find line starting with "Position of each atom"'
    natoms = len(types)

    while not line.startswith('Space Group '):
        line = f.readline()
        assert line, 'Cannot find line starting with "Space Group "'
    tokens = line.split()
    d['spacegroup'] = int(tokens[2])
    #d['symbol_nonconventional'] = tokens[3]
    d['symbol'] = tokens[4]

    d['origin'] = np.array([float(v) for v in f.readline().split()[2:]])
    f.readline()
    d['abc'] = np.array([[float(v) for v in f.readline().split()] 
                         for i in range(3)]).T

    while not line.startswith('Wyckoff position'):
        line = f.readline()
        assert line, 'Cannot find line starting with "Wyckoff position"'
    d['wyckoff'] = []
    d['tags'] = -np.ones(natoms, dtype=int)
    tag = 0
    while not line.startswith('---'):
        tokens = line.split()
        line = f.readline()
        assert line, 'Cannot find line starting with "Wyckoff position"'
        d['wyckoff'].append(tokens[2].rstrip(','))
        while not line.startswith('Wyckoff position') and not line.startswith('---'):
            tokens = line.split()
            n = int(tokens[0])
            d['tags'][n-1] = tag
            line = f.readline()
        tag += 1

    return d



def findsym(atoms, tol=1e-3, centering='P', types=None, isodata_dir=None):
    """Returns a dict describing the symmetry of *atoms*.

    Arguments
    ---------
    atoms: Atoms instance
        Atoms instance to find space group of.
    tol: float
        Accuracy to which dimensions of the unit cell and positions of
        atoms are known. Units in Angstrom.
    centering: 'P' | 'I' | 'F' | 'A' | 'B' | 'C' | 'R'
        Known centering: P (no known centering), I (body-centered), F
        (face-centered), A,B,C (base centered), R (rhombohedral
        centered with coordinates of centered points at (2/3,1/3,1/3)
        and (1/3,2/3,2/3)).
    types: None | sequence of integers
        Sequence of arbitrary positive integers identifying different
        atomic sites, so that a symmetry operation that takes one atom
        into another with different type would be forbidden.

    Returned dict items
    -------------------
    abc: 3x3 float array
        The vectors a, b, c defining the cell in scaled coordinates.
    cellpar: 6 floats
        Cell parameters a, b, c, alpha, beta, gamma with lengths in
        Angstrom and angles in degree.
    origin: 3 floats
        Origin of the space group with respect to the origin in the
        input data. Coordinates are dimensionless, given in terms of
        the lattice parameters of the unit cell in the input.
    spacegroup: int
        Space group number from the International Tables of
        Crystallography.
    symbol: str
        Hermann-Mauguin symbol (no spaces).
    tags: int array
        Array of site numbers for each atom.  Only atoms within the
        first conventional unit cell are tagged, the rest have -1 as
        tag.
    wyckoff: list
        List of wyckoff symbols for each site.
    """
    output = run(atoms, tol, centering, types, isodata_dir)
    d = parse(output)
    return d


def unique(atoms, tol=1e-3, centering='P', types=None, isodata_dir=None):
    """Returns an Atoms object containing only one atom from each unique site.
    """
    d = findsym(atoms, tol=tol, centering=centering, types=types, 
                isodata_dir=isodata_dir)
    mask = np.concatenate(([True], np.diff(d['tags']) != 0)) * (d['tags'] >= 0)
    at = atoms[mask]
    a, b, c, alpha, beta, gamma = d['cellpar']
    A, B, C = d['abc']
    A *= a
    B *= b
    C *= c
    from numpy.linalg import norm
    from numpy import cos, pi
    assert abs(np.dot(A, B) - (norm(A) * norm(B) * cos(gamma * pi / 180.))) < 1e-5
    assert abs(np.dot(A, C) - (norm(A) * norm(C) * cos(beta * pi / 180.))) < 1e-5
    assert abs(np.dot(B, C) - (norm(B) * norm(C) * cos(alpha * pi / 180.))) < 1e-5
    at.cell = np.array([A, B, C])
    for k in 'origin', 'spacegroup', 'wyckoff':
        at.info[k] = d[k]
    at.info['unit_cell'] = 'unique'
    scaled = at.get_scaled_positions()
    at.set_scaled_positions(scaled)
    return at


if __name__ == '__main__':
    import doctest
    print 'doctest: ', doctest.testmod()
