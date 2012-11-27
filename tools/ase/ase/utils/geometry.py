# Copyright (C) 2010, Jesper Friis
# (see accompanying license files for details).

"""Utility tools for convenient creation of slabs and interfaces of
different orientations."""

import numpy as np


def get_layers(atoms, miller, tolerance=0.001):
    """Returns two arrays describing which layer each atom belongs
    to and the distance between the layers and origo. 

    Parameters:

    miller: 3 integers
        The Miller indices of the planes. Actually, any direction
        in reciprocal space works, so if a and b are two float
        vectors spanning an atomic plane, you can get all layers
        parallel to this with miller=np.cross(a,b).
    tolerance: float
        The maximum distance in Angstrom along the plane normal for
        counting two atoms as belonging to the same plane.

    Returns:

    tags: array of integres
        Array of layer indices for each atom.
    levels: array of floats
        Array of distances in Angstrom from each layer to origo.

    Example:

    >>> import numpy as np
    >>> from ase.lattice.spacegroup import crystal
    >>> atoms = crystal('Al', [(0,0,0)], spacegroup=225, cellpar=4.05)
    >>> np.round(atoms.positions, decimals=5)
    array([[ 0.   ,  0.   ,  0.   ],
           [ 0.   ,  2.025,  2.025],
           [ 2.025,  0.   ,  2.025],
           [ 2.025,  2.025,  0.   ]])
    >>> get_layers(atoms, (0,0,1))
    (array([0, 1, 1, 0]), array([ 0.   ,  2.025]))
    """
    miller = np.asarray(miller)

    metric = np.dot(atoms.cell, atoms.cell.T)
    c = np.linalg.solve(metric.T, miller.T).T
    miller_norm = np.sqrt(np.dot(c, miller))
    d = np.dot(atoms.get_scaled_positions(), miller)/miller_norm

    keys = np.argsort(d)
    ikeys = np.argsort(keys)
    mask = np.concatenate(([True], np.diff(d[keys]) > tolerance))
    tags = np.cumsum(mask)[ikeys]
    if tags.min() == 1:
        tags -= 1

    levels = d[keys][mask]
    return tags, levels




def cut(atoms, a=(1, 0, 0), b=(0, 1, 0), c=None, clength=None, 
        origo=(0, 0, 0), nlayers=None, extend=1.0, tolerance=0.01, 
        maxatoms=None):
    """Cuts out a cell defined by *a*, *b*, *c* and *origo* from a
    sufficiently repeated copy of *atoms*.

    Typically, this function is used to create slabs of different
    sizes and orientations. The vectors *a*, *b* and *c* are in scaled
    coordinates and defines the returned cell and should normally be
    integer-valued in order to end up with a periodic
    structure. However, for systems with sub-translations, like fcc,
    integer multiples of 1/2 or 1/3 might also make sence for some
    directions (and will be treated correctly).

    Parameters:
    
    atoms: Atoms instance
        This should correspond to a repeatable unit cell.
    a: int | 3 floats
        The a-vector in scaled coordinates of the cell to cut out. If
        integer, the a-vector will be the scaled vector from *origo* to the
        atom with index *a*.
    b: int | 3 floats
        The b-vector in scaled coordinates of the cell to cut out. If
        integer, the b-vector will be the scaled vector from *origo* to the
        atom with index *b*.
    c: None | int | 3 floats
        The c-vector in scaled coordinates of the cell to cut out. 
        if integer, the c-vector will be the scaled vector from *origo* to 
        the atom with index *c*. 
        If *None* it will be along cross(a, b) converted to real space
        and normalised with the cube root of the volume. Note that this
        in general is not perpendicular to a and b for non-cubic
        systems. For cubic systems however, this is redused to 
        c = cross(a, b).
    clength: None | float
        If not None, the length of the c-vector will be fixed to
        *clength* Angstroms. Should not be used together with
        *nlayers*.
    origo: int | 3 floats
        Position of origo of the new cell in scaled coordinates. If
        integer, the position of the atom with index *origo* is used.
    nlayers: None | int
        If *nlayers* is not *None*, the returned cell will have
        *nlayers* atomic layers in the c-direction.
    extend: 1 or 3 floats
        The *extend* argument scales the effective cell in which atoms
        will be included. It must either be three floats or a single
        float scaling all 3 directions.  By setting to a value just
        above one, e.g. 1.05, it is possible to all the corner and
        edge atoms in the returned cell.  This will of cause make the
        returned cell non-repeatable, but is very usefull for
        visualisation.
    tolerance: float
        Determines what is defined as a plane.  All atoms within
        *tolerance* Angstroms from a given plane will be considered to
        belong to that plane.
    maxatoms: None | int
        This option is used to auto-tune *tolerance* when *nlayers* is
        given for high zone axis systems.  For high zone axis one
        needs to reduce *tolerance* in order to distinguise the atomic
        planes, resulting in the more atoms will be added and
        eventually MemoryError.  A too small *tolerance*, on the other
        hand, might result in inproper splitting of atomic planes and
        that too few layers are returned.  If *maxatoms* is not None,
        *tolerance* will automatically be gradually reduced until
        *nlayers* atomic layers is obtained, when the number of atoms
        exceeds *maxatoms*.

    Example:

    >>> import ase
    >>> from ase.lattice.spacegroup import crystal
    >>>
    # Create an aluminium (111) slab with three layers
    #
    # First an unit cell of Al
    >>> a = 4.05
    >>> aluminium = crystal('Al', [(0,0,0)], spacegroup=225,
    ...                     cellpar=[a, a, a, 90, 90, 90])
    >>>
    # Then cut out the slab
    >>> al111 = cut(aluminium, (1,-1,0), (0,1,-1), nlayers=3)
    >>>
    # Visualisation of the skutterudite unit cell 
    #
    # Again, create a skutterudite unit cell
    >>> a = 9.04
    >>> skutterudite = crystal(
    ...     ('Co', 'Sb'), 
    ...     basis=[(0.25,0.25,0.25), (0.0, 0.335, 0.158)], 
    ...     spacegroup=204, 
    ...     cellpar=[a, a, a, 90, 90, 90])
    >>>
    # Then use *origo* to put 'Co' at the corners and *extend* to
    # include all corner and edge atoms.
    >>> s = cut(skutterudite, origo=(0.25, 0.25, 0.25), extend=1.01)
    >>> ase.view(s)  # doctest: +SKIP
    """
    atoms = atoms.copy()
    cell = atoms.cell

    if isinstance(origo, int):
        origo = atoms.get_scaled_positions()[origo]
    origo = np.array(origo, dtype=float)

    scaled = (atoms.get_scaled_positions() - origo)%1.0
    scaled %= 1.0 # needed to ensure that all numbers are *less* than one
    atoms.set_scaled_positions(scaled)

    if isinstance(a, int):
        a = scaled[a] - origo
    if isinstance(b, int):
        b = scaled[b] - origo
    if isinstance(c, int):
        c = scaled[c] - origo

    a = np.array(a, dtype=float)
    b = np.array(b, dtype=float)
    if c is None:
        metric = np.dot(cell, cell.T)
        vol = np.sqrt(np.linalg.det(metric))
        h = np.cross(a, b) 
        H = np.linalg.solve(metric.T, h.T)
        c = vol*H/vol**(1./3.)
    c = np.array(c, dtype=float)

    if nlayers:
        # Recursive increase the length of c until we have at least
        # *nlayers* atomic layers parallell to the a-b plane
        while True:
            at = cut(atoms, a, b, c, origo=origo, extend=extend, 
                        tolerance=tolerance)
            scaled = at.get_scaled_positions()
            d = scaled[:,2]
            keys = np.argsort(d)
            ikeys = np.argsort(keys)
            tol = tolerance
            while True:
                mask = np.concatenate(([True], np.diff(d[keys]) > tol))
                tags = np.cumsum(mask)[ikeys] - 1
                levels = d[keys][mask]
                if (maxatoms is None or len(at) < maxatoms or 
                    len(levels) > nlayers): 
                    break
                tol *= 0.9
            if len(levels) > nlayers: 
                break
            c *= 2

        at.cell[2] *= levels[nlayers]
        return at[tags < nlayers]

    newcell = np.dot(np.array([a, b, c]), cell)
    if nlayers is None and clength is not None:
        newcell[2,:] *= clength/np.linalg.norm(newcell[2])

    # Create a new atoms object, repeated and translated such that
    # it completely covers the new cell
    scorners_newcell = np.array([[0., 0., 0.], [0., 0., 1.], 
                                 [0., 1., 0.], [0., 1., 1.], 
                                 [1., 0., 0.], [1., 0., 1.], 
                                 [1., 1., 0.], [1., 1., 1.]])
    corners = np.dot(scorners_newcell, newcell*extend)
    scorners = np.linalg.solve(cell.T, corners.T).T
    rep = np.ceil(scorners.ptp(axis=0)).astype('int') + 1
    trans = np.dot(np.floor(scorners.min(axis=0)), cell)
    atoms = atoms.repeat(rep)
    atoms.translate(trans)
    atoms.set_cell(newcell)

    # Mask out atoms outside new cell
    stol = 0.1*tolerance  # scaled tolerance, XXX
    maskcell = atoms.cell*extend
    sp = np.linalg.solve(maskcell.T, (atoms.positions).T).T
    mask = np.all(np.logical_and(-stol <= sp, sp < 1-stol), axis=1)
    atoms = atoms[mask]
    return atoms


class IncompatibleCellError(ValueError):
    """Exception raised if stacking fails due to incompatible cells
    between *atoms1* and *atoms2*."""
    pass


def stack(atoms1, atoms2, axis=2, cell=None, fix=0.5,  
          maxstrain=0.5, distance=None, reorder=False):
    """Return a new Atoms instance with *atoms2* stacked on top of
    *atoms1* along the given axis. Periodicity in all directions is
    ensured.

    The size of the final cell is determined by *cell*, except
    that the length alongh *axis* will be the sum of
    *atoms1.cell[axis]* and *atoms2.cell[axis]*. If *cell* is None,
    it will be interpolated between *atoms1* and *atoms2*, where
    *fix* determines their relative weight. Hence, if *fix* equals
    zero, the final cell will be determined purely from *atoms1* and
    if *fix* equals one, it will be determined purely from
    *atoms2*.

    An ase.geometry.IncompatibleCellError exception is raised if the
    cells of *atoms1* and *atoms2* are incopatible, e.g. if the far
    corner of the unit cell of either *atoms1* or *atoms2* is
    displaced more than *maxstrain*. Setting *maxstrain* to None,
    disable this check.

    If *distance* is not None, the size of the final cell, along the
    direction perpendicular to the interface, will be adjusted such
    that the distance between the closest atoms in *atoms1* and
    *atoms2* will be equal to *distance*. This option uses
    scipy.optimize.fmin() and hence require scipy to be installed.

    If *reorder* is True, then the atoms will be reordred such that
    all atoms with the same symbol will follow sequensially after each
    other, eg: 'Al2MnAl10Fe' -> 'Al12FeMn'.    

    Example:

    >>> import ase
    >>> from ase.lattice.spacegroup import crystal
    >>>
    # Create an Ag(110)-Si(110) interface with three atomic layers
    # on each side. 
    >>> a_ag = 4.09
    >>> ag = crystal(['Ag'], basis=[(0,0,0)], spacegroup=225, 
    ...              cellpar=[a_ag, a_ag, a_ag, 90., 90., 90.])
    >>> ag110 = cut(ag, (0, 0, 3), (-1.5, 1.5, 0), nlayers=3)
    >>>
    >>> a_si = 5.43
    >>> si = crystal(['Si'], basis=[(0,0,0)], spacegroup=227, 
    ...              cellpar=[a_si, a_si, a_si, 90., 90., 90.])
    >>> si110 = cut(si, (0, 0, 2), (-1, 1, 0), nlayers=3)
    >>>
    >>> interface = stack(ag110, si110, maxstrain=1)
    >>> ase.view(interface)  # doctest: +SKIP
    >>>
    # Once more, this time adjusted such that the distance between
    # the closest Ag and Si atoms will be 2.3 Angstrom (requires scipy).
    >>> interface2 = stack(ag110, si110, 
    ...                    maxstrain=1, distance=2.3)   # doctest:+ELLIPSIS
    Optimization terminated successfully.
        ...
    >>> ase.view(interface2)  # doctest: +SKIP
    """
    atoms1 = atoms1.copy()
    atoms2 = atoms2.copy()

    if (np.sign(np.linalg.det(atoms1.cell)) != 
        np.sign(np.linalg.det(atoms2.cell))):
        raise IncompatibleCellError('*atoms1* amd *atoms2* must both either '
                                    'have a lefthanded or a righanded cell.')

    c1 = np.linalg.norm(atoms1.cell[axis])
    c2 = np.linalg.norm(atoms2.cell[axis])
    if cell is None:
        cell1 = atoms1.cell.copy()
        cell2 = atoms2.cell.copy()
        cell1[axis] /= c1
        cell2[axis] /= c2
        cell = cell1 + fix*(cell2 - cell1)
    cell[axis] /= np.linalg.norm(cell[axis])
    cell1 = cell.copy()
    cell2 = cell.copy()
    cell1[axis] *= c1
    cell2[axis] *= c2

    if maxstrain:
        strain1 = np.sqrt(((cell1 - atoms1.cell).sum(axis=0)**2).sum())
        strain2 = np.sqrt(((cell2 - atoms2.cell).sum(axis=0)**2).sum())
        if strain1 > maxstrain or strain2 > maxstrain:
            raise IncompatibleCellError(
                '*maxstrain* exceeded. *atoms1* strained %f and '
                '*atoms2* strained %f.'%(strain1, strain2))

    atoms1.set_cell(cell1, scale_atoms=True)
    atoms2.set_cell(cell2, scale_atoms=True)

    if distance is not None:
        from scipy.optimize import fmin
        def mindist(pos1, pos2):
            n1 = len(pos1)
            n2 = len(pos2)
            idx1 = np.arange(n1).repeat(n2)
            idx2 = np.tile(np.arange(n2), n1)
            return np.sqrt(((pos1[idx1] - pos2[idx2])**2).sum(axis=1).min())
        def func(x):
            t1, t2, h1, h2 = x[0:3], x[3:6], x[6], x[7]
            pos1 = atoms1.positions + t1
            pos2 = atoms2.positions + t2
            d1 = mindist(pos1, pos2 + (h1 + 1.0)*atoms1.cell[axis])
            d2 = mindist(pos2, pos1 + (h2 + 1.0)*atoms2.cell[axis])
            return (d1 - distance)**2 + (d2 - distance)**2
        atoms1.center()
        atoms2.center()
        x0 = np.zeros((8,))
        x = fmin(func, x0)
        t1, t2, h1, h2 = x[0:3], x[3:6], x[6], x[7]
        atoms1.translate(t1)
        atoms2.translate(t2)
        atoms1.cell[axis] *= 1.0 + h1
        atoms2.cell[axis] *= 1.0 + h2

    atoms2.translate(atoms1.cell[axis])
    atoms1.cell[axis] += atoms2.cell[axis]
    atoms1.extend(atoms2)

    if reorder:
        atoms1 = sort(atoms1)

    return atoms1


def sort(atoms, tags=None):
    """Return a new Atoms object with sorted atomic order. The default
    is to order according to chemical symbols, but if *tags* is not
    None, it will be used instead. A stable sorting algorithm is used.

    Example:
    >>> import ase
    >>> from ase.lattice.spacegroup import crystal
    >>>
    # Two unit cells of NaCl
    >>> a = 5.64
    >>> nacl = crystal(['Na', 'Cl'], [(0, 0, 0), (0.5, 0.5, 0.5)], 
    ... spacegroup=225, cellpar=[a, a, a, 90, 90, 90]).repeat((2, 1, 1))
    >>> nacl.get_chemical_symbols()
    ['Na', 'Na', 'Na', 'Na', 'Cl', 'Cl', 'Cl', 'Cl', 'Na', 'Na', 'Na', 'Na', 'Cl', 'Cl', 'Cl', 'Cl']
    >>> nacl_sorted = sort(nacl)
    >>> nacl_sorted.get_chemical_symbols()
    ['Cl', 'Cl', 'Cl', 'Cl', 'Cl', 'Cl', 'Cl', 'Cl', 'Na', 'Na', 'Na', 'Na', 'Na', 'Na', 'Na', 'Na']
    >>> np.all(nacl_sorted.cell == nacl.cell)
    True
    """
    if tags is None:
        tags = atoms.get_chemical_symbols()
    else:
        tags = list(tags)
    deco = [(tag, i) for i, tag in enumerate(tags)]
    deco.sort()
    indices = [i for tag, i in deco]
    return atoms[indices]


def rotation_matrix(a1, a2, b1, b2):
    """Returns a rotation matrix that rotates the vectors *a1* in the
    direction of *a2* and *b1* in the direction of *b2*.
    
    In the case that the angle between *a2* and *b2* is not the same
    as between *a1* and *b1*, a proper rotation matrix will anyway be
    constructed by first rotate *b2* in the *b1*, *b2* plane.
    """
    from numpy.linalg import norm, det
    a1 = np.asarray(a1, dtype=float)/norm(a1)
    b1 = np.asarray(b1, dtype=float)/norm(b1)
    c1 = np.cross(a1, b1)
    c1 /= norm(c1)      # clean out rounding errors...

    a2 = np.asarray(a2, dtype=float)/norm(a2)
    b2 = np.asarray(b2, dtype=float)/norm(b2)
    c2 = np.cross(a2, b2)
    c2 /= norm(c2)      # clean out rounding errors...

    # Calculate rotated *b2*
    theta = np.arccos(np.dot(a2, b2)) - np.arccos(np.dot(a1, b1))
    b3 = np.sin(theta)*a2 + np.cos(theta)*b2
    b3 /= norm(b3)      # clean out rounding errors...
    
    A1 = np.array([a1, b1, c1])
    A2 = np.array([a2, b3, c2])
    R = np.linalg.solve(A1, A2).T
    return R


def rotate(atoms, a1, a2, b1, b2, rotate_cell=True, center=(0, 0, 0)):
    """Rotate *atoms*, such that *a1* will be rotated in the direction
    of *a2* and *b1* in the direction of *b2*.  The point at *center*
    is fixed.  Use *center='COM'* to fix the center of mass.  If
    *rotate_cell* is true, the cell will be rotated together with the
    atoms.

    Note that the 000-corner of the cell is by definition fixed at
    origo.  Hence, setting *center* to something other than (0, 0, 0)
    will rotate the atoms out of the cell, even if *rotate_cell* is
    True.
    """
    if isinstance(center, str) and center.lower() == 'com':
        center = atoms.get_center_of_mass()
    
    R = rotation_matrix(a1, a2, b1, b2)
    atoms.positions[:] = np.dot(atoms.positions - center, R.T) + center

    if rotate_cell:
        atoms.cell[:] = np.dot(atoms.cell, R.T)

#-----------------------------------------------------------------
# Self test
if __name__ == '__main__':
    import doctest
    print 'doctest: ', doctest.testmod()
