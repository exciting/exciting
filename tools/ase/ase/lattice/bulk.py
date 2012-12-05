from math import sqrt

from ase.atoms import Atoms, string2symbols
from ase.data import reference_states, atomic_numbers, chemical_symbols


def bulk(name, crystalstructure=None, a=None, c=None, covera=None,
         orthorhombic=False, cubic=False):
    """Creating bulk systems.

    Crystal structure and lattice constant(s) will be guessed if not
    provided.

    name: str
        Chemical symbol or symbols as in 'MgO' or 'NaCl'.
    crystalstructure: str
        Must be one of sc, fcc, bcc, hcp, diamond, zincblende,
        rocksalt, cesiumchloride, or fluorite.
    a: float
        Lattice constant.
    c: float
        Lattice constant.
    covera: float
        c/a raitio used for hcp.  Use sqrt(8/3.0) for ideal ratio.
    orthorhombic: bool
        Construct orthorhombic unit cell instead of primitive cell
        which is the default.
    cubic: bool
        Construct cubic unit cell if possible.
    """

    if a is not None:
        a = float(a)
    if c is not None:
        c = float(c)

    if covera is not None and c is not None:
        raise ValueError("Don't specify both c and c/a!")

    if name in chemical_symbols:
        Z = atomic_numbers[name]
        ref = reference_states[Z]
        if ref is not None:
            xref = ref['symmetry']
        else:
            xref = None

    if crystalstructure is None:
        crystalstructure = xref

    if a is None:
        assert xref == crystalstructure
        a = ref['a']

    if crystalstructure == 'hcp':
        cubic = False
        if c is not None:
            covera = c / a
        elif covera is None:
            if xref == 'hcp':
                covera = ref['c/a']
            else:
                covera = sqrt(8.0 / 3.0)

    if orthorhombic and crystalstructure != 'sc':
        return _orthorhombic_bulk(name, crystalstructure, a, covera)

    if cubic and crystalstructure in ['bcc', 'cesiumchloride']:
        return _orthorhombic_bulk(name, crystalstructure, a, covera)

    if cubic and crystalstructure != 'sc':
        return _cubic_bulk(name, crystalstructure, a)

    if crystalstructure == 'sc':
        atoms = Atoms(name, cell=(a, a, a), pbc=True)
    elif crystalstructure == 'fcc':
        b = a / 2
        atoms = Atoms(name, cell=[(0, b, b), (b, 0, b), (b, b, 0)], pbc=True)
    elif crystalstructure == 'bcc':
        b = a / 2
        atoms = Atoms(name, cell=[(-b, b, b), (b, -b, b), (b, b, -b)],
                      pbc=True)
    elif crystalstructure == 'hcp':
        atoms = Atoms(2 * name,
                      scaled_positions=[(0, 0, 0),
                                        (1.0 / 3.0, 2.0 / 3.0, 0.5)],
                      cell=[(a, 0, 0),
                            (-a / 2, a * sqrt(3) / 2, 0),
                            (0, 0, covera * a)],
                      pbc=True)
    elif crystalstructure == 'diamond':
        atoms = bulk(2 * name, 'zincblende', a)
    elif crystalstructure == 'zincblende':
        s1, s2 = string2symbols(name)
        atoms = bulk(s1, 'fcc', a) + bulk(s2, 'fcc', a)
        atoms.positions[1] += a / 4
    elif crystalstructure == 'rocksalt':
        s1, s2 = string2symbols(name)
        atoms = bulk(s1, 'fcc', a) + bulk(s2, 'fcc', a)
        atoms.positions[1, 0] += a / 2
    elif crystalstructure == 'cesiumchloride':
        s1, s2 = string2symbols(name)
        atoms = bulk(s1, 'sc', a) + bulk(s2, 'sc', a)
        atoms.positions[1, :] += a / 2
    elif crystalstructure == 'fluorite':
        s1, s2, s3 = string2symbols(name)
        atoms = bulk(s1, 'fcc', a) + bulk(s2, 'fcc', a) + bulk(s3, 'fcc', a)
        atoms.positions[1, :] += a / 4
        atoms.positions[2, :] += a * 3. / 4
    else:
        raise ValueError('Unknown crystal structure: ' + crystalstructure)

    return atoms


def _orthorhombic_bulk(name, crystalstructure, a, covera=None):
    if crystalstructure == 'fcc':
        b = a / sqrt(2)
        atoms = Atoms(2 * name, cell=(b, b, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)])
    elif crystalstructure == 'bcc':
        atoms = Atoms(2 * name, cell=(a, a, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)])
    elif crystalstructure == 'hcp':
        atoms = Atoms(4 * name,
                      cell=(a, a * sqrt(3), covera * a),
                      scaled_positions=[(0, 0, 0),
                                        (0.5, 0.5, 0),
                                        (0.5, 1.0 / 6.0, 0.5),
                                        (0, 2.0 / 3.0, 0.5)],
                      pbc=True)
    elif crystalstructure == 'diamond':
        atoms = _orthorhombic_bulk(2 * name, 'zincblende', a)
    elif crystalstructure == 'zincblende':
        s1, s2 = string2symbols(name)
        b = a / sqrt(2)
        atoms = Atoms(2 * name, cell=(b, b, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0, 0.25),
                                        (0.5, 0.5, 0.5), (0, 0.5, 0.75)])
    elif crystalstructure == 'rocksalt':
        s1, s2 = string2symbols(name)
        b = a / sqrt(2)
        atoms = Atoms(2 * name, cell=(b, b, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0.5, 0),
                                        (0.5, 0.5, 0.5), (0, 0, 0.5)])
    elif crystalstructure == 'cesiumchloride':
        atoms = Atoms(name, cell=(a, a, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5)])
    else:
        raise RuntimeError

    return atoms


def _cubic_bulk(name, crystalstructure, a):
    if crystalstructure == 'fcc':
        atoms = Atoms(4 * name, cell=(a, a, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0, 0.5, 0.5),
                                        (0.5, 0, 0.5), (0.5, 0.5, 0)])
    elif crystalstructure == 'diamond':
        atoms = _cubic_bulk(2 * name, 'zincblende', a)
    elif crystalstructure == 'zincblende':
        atoms = Atoms(4 * name, cell=(a, a, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.25, 0.25, 0.25),
                                        (0, 0.5, 0.5), (0.25, 0.75, 0.75),
                                        (0.5, 0, 0.5), (0.75, 0.25, 0.75),
                                        (0.5, 0.5, 0), (0.75, 0.75, 0.25)])
    elif crystalstructure == 'rocksalt':
        atoms = Atoms(4 * name, cell=(a, a, a), pbc=True,
                      scaled_positions=[(0, 0, 0), (0.5, 0, 0),
                                        (0, 0.5, 0.5), (0.5, 0.5, 0.5),
                                        (0.5, 0, 0.5), (0, 0, 0.5),
                                        (0.5, 0.5, 0), (0, 0.5, 0)])
    else:
        raise RuntimeError

    return atoms
