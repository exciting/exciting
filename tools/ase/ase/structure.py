"""Atomic structure.

This mudule contains helper functions for setting up nanotubes and
graphene nanoribbons."""

import warnings
from math import sqrt

import numpy as np

from ase.atoms import Atoms, string2symbols
from ase.data import covalent_radii
from ase.utils import gcd


def nanotube(n, m, length=1, bond=1.42, symbol='C', verbose=False):
    if n < m:
        m, n = n, m
        sign = -1
    else:
        sign = 1

    nk = 6000
    sq3 = sqrt(3.0)
    a = sq3 * bond
    l2 = n * n + m * m + n * m
    l = sqrt(l2)
    dt = a * l / np.pi

    nd = gcd(n ,m)
    if (n - m) % (3 * nd ) == 0:
        ndr = 3 * nd
    else:
        ndr = nd

    nr = (2 * m + n) / ndr
    ns = -(2 * n + m) / ndr
    nt2 = 3 * l2 / ndr / ndr
    nt = np.floor(sqrt(nt2))
    nn = 2 * l2 / ndr

    ichk = 0
    if nr == 0:
        n60 = 1
    else:
        n60 = nr * 4

    absn = abs(n60)
    nnp = []
    nnq = []
    for i in range(-absn, absn + 1):
        for j in range(-absn, absn + 1):
            j2 = nr * j - ns * i
            if j2 == 1:
                j1 = m * i - n * j
                if j1 > 0 and j1 < nn:
                    ichk += 1
                    nnp.append(i)
                    nnq.append(j)

    if ichk == 0:
        raise RuntimeError('not found p, q strange!!')
    if ichk >= 2:
        raise RuntimeError('more than 1 pair p, q strange!!')

    nnnp = nnp[0]
    nnnq = nnq[0]

    if verbose:
        print 'the symmetry vector is', nnnp, nnnq

    lp = nnnp * nnnp + nnnq * nnnq + nnnp * nnnq
    r = a * sqrt(lp)
    c = a * l
    t = sq3 * c / ndr

    if 2 * nn > nk:
        raise RuntimeError('parameter nk is too small!')

    rs = c / (2.0 * np.pi)

    if verbose:
        print 'radius=', rs, t

    q1 = np.arctan((sq3 * m) / (2 * n + m))
    q2 = np.arctan((sq3 * nnnq) / (2 * nnnp + nnnq))
    q3 = q1 - q2

    q4 = 2.0 * np.pi / nn
    q5 = bond * np.cos((np.pi / 6.0) - q1) / c * 2.0 * np.pi

    h1 = abs(t) / abs(np.sin(q3))
    h2 = bond * np.sin((np.pi / 6.0) - q1)

    ii = 0
    x, y, z = [], [], []
    for i in range(nn):
        x1, y1, z1 = 0, 0, 0

        k = np.floor(i * abs(r) / h1)
        x1 = rs * np.cos(i * q4)
        y1 = rs * np.sin(i * q4)
        z1 = (i * abs(r) - k * h1) * np.sin(q3)
        kk2 = abs(np.floor((z1 + 0.0001) / t))
        if z1 >= t - 0.0001:
            z1 -= t * kk2
        elif z1 < 0:
            z1 += t * kk2
        ii += 1

        x.append(x1)
        y.append(y1)
        z.append(z1)
        z3 = (i * abs(r) - k * h1) * np.sin(q3) - h2
        ii += 1

        if z3 >= 0 and z3 < t:
            x2 = rs * np.cos(i * q4 + q5)
            y2 = rs * np.sin(i * q4 + q5)
            z2 = (i * abs(r) - k * h1) * np.sin(q3) - h2
            x.append(x2)
            y.append(y2)
            z.append(z2)
        else:
            x2 = rs * np.cos(i * q4 + q5)
            y2 = rs * np.sin(i * q4 + q5)
            z2 = (i * abs(r) - (k + 1) * h1) * np.sin(q3) - h2
            kk = abs(np.floor(z2 / t))
            if z2 >= t - 0.0001:
                z2 -= t * kk
            elif z2 < 0:
                z2 += t * kk
            x.append(x2)
            y.append(y2)
            z.append(z2)

    ntotal = 2 * nn
    X = []
    for i in range(ntotal):
        X.append([x[i], y[i], sign * z[i]])

    if length > 1:
        xx = X[:]
        for mnp in range(2, length + 1):
            for i in range(len(xx)):
                X.append(xx[i][:2] + [xx[i][2] + (mnp - 1) * t])

    TransVec = t
    NumAtom = ntotal * length
    Diameter = rs * 2
    ChiralAngle = np.arctan((sq3 * n) / (2 * m + n)) / (np.pi * 180)

    cell = [Diameter * 2, Diameter * 2, length * t]
    atoms = Atoms(symbol + str(NumAtom), positions=X, cell=cell,
                  pbc=[False, False, True])
    atoms.center()
    if verbose:
        print 'translation vector =', TransVec
        print 'diameter = ', Diameter
        print 'chiral angle = ', ChiralAngle
    return atoms

def graphene_nanoribbon(n, m, type='zigzag', saturated=False, C_H=1.09,
                        C_C=1.42, vacuum=2.5, magnetic=None, initial_mag=1.12,
                        sheet=False, main_element='C', saturate_element='H',
                        vacc=None):
    """Create a graphene nanoribbon.

    Creates a graphene nanoribbon in the x-z plane, with the nanoribbon
    running along the z axis.

    Parameters:

    n: The width of the nanoribbon

    m: The length of the nanoribbon.

    type ('zigzag'): The orientation of the ribbon.  Must be either 'zigzag'
    or 'armchair'.

    saturated (Falsi):  If true, hydrogen atoms are placed along the edge.

    C_H: Carbon-hydrogen bond length.  Default: 1.09 Angstrom

    C_C: Carbon-carbon bond length.  Default: 1.42 Angstrom.

    vacuum:  Amount of vacuum added to both sides.  Default 2.5 Angstrom.

    magnetic:  Make the edges magnetic.

    initial_mag: Magnitude of magnetic moment if magnetic=True.

    sheet:  If true, make an infinite sheet instead of a ribbon.
    """
    #This function creates the coordinates for a graphene nanoribbon,
    #n is width, m is length

    if vacc is not None:
        warnings.warn('Use vacuum=%f' % (0.5 * vacc))
        vacuum = 0.5 * vacc

    assert vacuum > 0
    b = sqrt(3) * C_C / 4
    arm_unit = Atoms(main_element+'4', pbc=(1,0,1),
                     cell = [4 * b,  2 * vacuum,  3 * C_C])
    arm_unit.positions = [[0, 0, 0],
                          [b * 2, 0, C_C / 2.],
                          [b * 2, 0, 3 * C_C / 2.],
                          [0, 0, 2 * C_C]]
    zz_unit = Atoms(main_element+'2', pbc=(1,0,1),
                    cell = [3 * C_C /2., 2 * vacuum, b * 4])
    zz_unit.positions = [[0, 0, 0],
                         [C_C / 2., 0, b * 2]]
    atoms = Atoms()
    tol = 1e-4
    if sheet:
        vacuum2 = 0.0
    else:
        vacuum2 = vacuum
    if type == 'zigzag':
        edge_index0 = np.arange(m) * 2 + 1
        edge_index1 = (n - 1) * m * 2 + np.arange(m) * 2
        if magnetic:
            mms = np.zeros(m * n * 2)
            for i in edge_index0:
                mms[i] = initial_mag
            for i in edge_index1:
                mms[i] = -initial_mag

        for i in range(n):
            layer = zz_unit.repeat((1, 1, m))
            layer.positions[:, 0] -= 3 * C_C / 2 * i
            if i % 2 == 1:
                layer.positions[:, 2] += 2 * b
                layer[-1].position[2] -= b * 4 * m
            atoms += layer
        if magnetic:
            atoms.set_initial_magnetic_moments(mms)
        if saturated:
            H_atoms0 = Atoms(saturate_element + str(m))
            H_atoms0.positions = atoms[edge_index0].positions
            H_atoms0.positions[:, 0] += C_H
            H_atoms1 = Atoms(saturate_element + str(m))
            H_atoms1.positions = atoms[edge_index1].positions
            H_atoms1.positions[:, 0] -= C_H
            atoms += H_atoms0 + H_atoms1
        atoms.cell = [n * 3 * C_C / 2 + 2 * vacuum2, 2 * vacuum, m * 4 * b]

    elif type == 'armchair':
        for i in range(n):
            layer = arm_unit.repeat((1, 1, m))
            layer.positions[:, 0] -= 4 * b * i
            atoms += layer
        atoms.cell = [b * 4 * n + 2 * vacuum2, 2 * vacuum, 3 * C_C * m]

    atoms.center()
    atoms.set_pbc([sheet, False, True])
    return atoms

def molecule(name, data=None, **kwargs):
    """Create formula base on data. If data is None assume G2 set.
    kwargs currently not used.  """
    if data is None:
        from ase.data.g2 import data
    if name not in data.keys():
        raise NotImplementedError('%s not in data.' % (name))
    args = data[name].copy()
    # accept only the following Atoms constructor arguments
    # XXX: should we accept all Atoms arguments?
    for k in args.keys():
        if k not in [
            'symbols', 'positions', 'numbers',
            'tags', 'masses',
            'magmoms', 'charges',
            'info',
            ]:
            args.pop(k)
    # kwargs overwrites data
    args.update(kwargs)
    return Atoms(**args)

def bulk(name, crystalstructure, a=None, c=None, covera=None,
         orthorhombic=False, cubic=False):
    """Helper function for creating bulk systems.

    name: str
        Chemical symbol or symbols as in 'MgO' or 'NaCl'.
    crystalstructure: str
        Must be one of sc, fcc, bcc, hcp, diamond, zincblende or
        rocksalt.
    a: float
        Lattice constant.
    c: float
        Lattice constant.
    covera: float
        c/a raitio used for hcp.  Defaults to ideal ratio.
    orthorhombic: bool
        Construct orthorhombic unit cell instead of primitive cell
        which is the default.
    cubic: bool
        Construct cubic unit cell.
    """

    #warnings.warn('This function is deprecated.  Use the ' +
    #              'ase.lattice.bulk.bulk() function instead.')

    if a is not None:
        a = float(a)
    if c is not None:
        c = float(c)

    if covera is not None and  c is not None:
        raise ValueError("Don't specify both c and c/a!")

    if covera is None and c is None:
        covera = sqrt(8.0 / 3.0)

    if a is None:
        a = estimate_lattice_constant(name, crystalstructure, covera)

    if covera is None and c is not None:
        covera = c / a

    x = crystalstructure.lower()

    if orthorhombic and x != 'sc':
        return _orthorhombic_bulk(name, x, a, covera)

    if cubic and x == 'bcc':
        return _orthorhombic_bulk(name, x, a, covera)

    if cubic and x != 'sc':
        return _cubic_bulk(name, x, a)

    if x == 'sc':
        atoms = Atoms(name, cell=(a, a, a), pbc=True)
    elif x == 'fcc':
        b = a / 2
        atoms = Atoms(name, cell=[(0, b, b), (b, 0, b), (b, b, 0)], pbc=True)
    elif x == 'bcc':
        b = a / 2
        atoms = Atoms(name, cell=[(-b, b, b), (b, -b, b), (b, b, -b)],
                      pbc=True)
    elif x == 'hcp':
        atoms = Atoms(2 * name,
                      scaled_positions=[(0, 0, 0),
                                        (1.0 / 3.0, 1.0 / 3.0, 0.5)],
                      cell=[(a, 0, 0),
                            (a / 2, a * sqrt(3) / 2, 0),
                            (0, 0, covera * a)],
                      pbc=True)
    elif x == 'diamond':
        atoms = bulk(2 * name, 'zincblende', a)
    elif x == 'zincblende':
        s1, s2 = string2symbols(name)
        atoms = bulk(s1, 'fcc', a) + bulk(s2, 'fcc', a)
        atoms.positions[1] += a / 4
    elif x == 'rocksalt':
        s1, s2 = string2symbols(name)
        atoms = bulk(s1, 'fcc', a) + bulk(s2, 'fcc', a)
        atoms.positions[1, 0] += a / 2
    else:
        raise ValueError('Unknown crystal structure: ' + crystalstructure)

    return atoms

def estimate_lattice_constant(name, crystalstructure, covera):
    atoms = bulk(name, crystalstructure, 1.0, covera)
    v0 = atoms.get_volume()
    v = 0.0
    for Z in atoms.get_atomic_numbers():
        r = covalent_radii[Z]
        v += 4 * np.pi / 3 * r**3 * 1.5
    return (v / v0)**(1.0 / 3)

def _orthorhombic_bulk(name, x, a, covera=None):
    if x == 'fcc':
        b = a / sqrt(2)
        atoms = Atoms(2 * name, cell=(b, b, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)])
    elif x == 'bcc':
        atoms = Atoms(2 * name, cell=(a, a, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)])
    elif x == 'hcp':
        atoms = Atoms(4 * name,
                      cell=(a, a * sqrt(3), covera * a),
                      scaled_positions=[(0, 0, 0),
                                        (0.5, 0.5, 0),
                                        (0.5, 1.0 / 6.0, 0.5),
                                        (0, 2.0 / 3.0, 0.5)],
                      pbc=True)
    elif x == 'diamond':
        atoms = _orthorhombic_bulk(2 * name, 'zincblende', a)
    elif x == 'zincblende':
        s1, s2 = string2symbols(name)
        b = a / sqrt(2)
        atoms = Atoms(2 * name, cell=(b, b, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0, 0.25),
                                        (0.5, 0.5, 0.5), (0, 0.5, 0.75)])
    elif x == 'rocksalt':
        s1, s2 = string2symbols(name)
        b = a / sqrt(2)
        atoms = Atoms(2 * name, cell=(b, b, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0.5, 0),
                                        (0.5, 0.5, 0.5), (0, 0, 0.5)])
    else:
        raise RuntimeError

    return atoms

def _cubic_bulk(name, x, a):
    if x == 'fcc':
        atoms = Atoms(4 * name, cell=(a, a, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0, 0.5, 0.5),
                                        (0.5, 0, 0.5), (0.5, 0.5, 0)])
    elif x == 'diamond':
        atoms = _cubic_bulk(2 * name, 'zincblende', a)
    elif x == 'zincblende':
        atoms = Atoms(4 * name, cell=(a, a, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.25, 0.25, 0.25),
                                        (0, 0.5, 0.5), (0.25, 0.75, 0.75),
                                        (0.5, 0, 0.5), (0.75, 0.25, 0.75),
                                        (0.5, 0.5, 0), (0.75, 0.75, 0.25)])
    elif x == 'rocksalt':
        atoms = Atoms(4 * name, cell=(a, a, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0, 0),
                                        (0, 0.5, 0.5), (0.5, 0.5, 0.5),
                                        (0.5, 0, 0.5), (0, 0, 0.5),
                                        (0.5, 0.5, 0), (0, 0.5, 0)])
    else:
        raise RuntimeError

    return atoms
