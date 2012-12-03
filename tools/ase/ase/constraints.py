from math import sqrt

import numpy as np

__all__ = ['FixCartesian', 'FixBondLength', 'FixedMode', 'FixConstraintSingle',
           'FixAtoms', 'UnitCellFilter', 'FixScaled', 'StrainFilter',
           'FixedPlane', 'Filter', 'FixConstraint', 'FixedLine',
           'FixBondLengths', 'FixInternals']


def slice2enlist(s):
    """Convert a slice object into a list of (new, old) tuples."""
    if isinstance(s, (list, tuple)):
        return enumerate(s)
    if s.step == None:
        step = 1
    else:
        step = s.step
    if s.start == None:
        start = 0
    else:
        start = s.start
    return enumerate(range(start, s.stop, step))


class FixConstraint:
    """Base class for classes that fix one or more atoms in some way."""

    def index_shuffle(self, ind):
        """Change the indices.

        When the ordering of the atoms in the Atoms object changes,
        this method can be called to shuffle the indices of the
        constraints.

        ind -- List or tuple of indices.

        """
        raise NotImplementedError

    def repeat(self, m, n):
        """ basic method to multiply by m, needs to know the length
        of the underlying atoms object for the assignment of
        multiplied constraints to work.
        """
        raise NotImplementedError


class FixConstraintSingle(FixConstraint):
    """Base class for classes that fix a single atom."""

    def index_shuffle(self, ind):
        """The atom index must be stored as self.a."""
        newa = -1   # Signal error
        for new, old in slice2enlist(ind):
            if old == self.a:
                newa = new
                break
        if newa == -1:
            raise IndexError('Constraint not part of slice')
        self.a = newa


class FixAtoms(FixConstraint):
    """Constraint object for fixing some chosen atoms."""
    def __init__(self, indices=None, mask=None):
        """Constrain chosen atoms.

        Parameters
        ----------
        indices : list of int
           Indices for those atoms that should be constrained.
        mask : list of bool
           One boolean per atom indicating if the atom should be
           constrained or not.

        Examples
        --------
        Fix all Copper atoms:

        >>> c = FixAtoms(mask=[s == 'Cu' for s in atoms.get_chemical_symbols()])
        >>> atoms.set_constraint(c)

        Fix all atoms with z-coordinate less than 1.0 Angstrom:

        >>> c = FixAtoms(mask=atoms.positions[:, 2] < 1.0)
        >>> atoms.set_constraint(c)
        """

        if indices is None and mask is None:
            raise ValueError('Use "indices" or "mask".')
        if indices is not None and mask is not None:
            raise ValueError('Use only one of "indices" and "mask".')

        if mask is not None:
            self.index = np.asarray(mask, bool)
        else:
            # Check for duplicates
            srt = np.sort(indices)
            for i in range(len(indices) - 1):
                if srt[i] == srt[i+1]:
                    raise ValueError(
                        'FixAtoms: The indices array contained duplicates. '
                        'Perhaps you wanted to specify a mask instead, but '
                        'forgot the mask= keyword.')
            self.index = np.asarray(indices, int)

        if self.index.ndim != 1:
            raise ValueError('Wrong argument to FixAtoms class!')

    def adjust_positions(self, old, new):
        new[self.index] = old[self.index]

    def adjust_forces(self, positions, forces):
        forces[self.index] = 0.0

    def index_shuffle(self, ind):
        # See docstring of superclass
        if self.index.dtype == bool:
            self.index = self.index[ind]
        else:
            index = []
            for new, old in slice2enlist(ind):
                if old in self.index:
                    index.append(new)
            if len(index) == 0:
                raise IndexError('All indices in FixAtoms not part of slice')
            self.index = np.asarray(index, int)

    def copy(self):
        if self.index.dtype == bool:
            return FixAtoms(mask=self.index.copy())
        else:
            return FixAtoms(indices=self.index.copy())

    def __repr__(self):
        if self.index.dtype == bool:
            return 'FixAtoms(mask=%s)' % ints2string(self.index.astype(int))
        return 'FixAtoms(indices=%s)' % ints2string(self.index)

    def repeat(self, m, n):
        i0 = 0
        l = len(self.index)
        natoms = 0
        if isinstance(m, int):
            m = (m, m, m)
        index_new = []
        for m2 in range(m[2]):
            for m1 in range(m[1]):
                for m0 in range(m[0]):
                    i1 = i0 + n
                    if self.index.dtype == bool:
                        index_new.extend(self.index)
                    else:
                        index_new += [i+natoms for i in self.index]
                    i0 = i1
                    natoms += n
        if self.index.dtype == bool:
            self.index = np.asarray(index_new, bool)
        else:
            self.index = np.asarray(index_new, int)
        return self

    def delete_atom(self, ind):
        """ Removes atom number ind from the index array, if present.
        Required for removing atoms with existing FixAtoms constraints.
        """
        if self.index.dtype == bool:
            self.index = np.delete(self.index, ind)
        else:
            if ind in self.index:
                i = list(self.index).index(ind)
                self.index = np.delete(self.index, i)
            for i in range(len(self.index)):
                if self.index[i] >= ind:
                    self.index[i] -= 1

def ints2string(x, threshold=10):
    """Convert ndarray of ints to string."""
    if len(x) <= threshold:
        return str(x.tolist())
    return str(x[:threshold].tolist())[:-1] + ', ...]'

class FixBondLengths(FixConstraint):
    def __init__(self, pairs, iterations=10):
        self.constraints = [FixBondLength(a1, a2)
                            for a1, a2 in pairs]
        self.iterations = iterations

    def adjust_positions(self, old, new):
        for i in range(self.iterations):
            for constraint in self.constraints:
                constraint.adjust_positions(old, new)

    def adjust_forces(self, positions, forces):
        for i in range(self.iterations):
            for constraint in self.constraints:
                constraint.adjust_forces(positions, forces)

    def copy(self):
        return FixBondLengths([constraint.indices
                               for constraint in self.constraints])

class FixBondLength(FixConstraint):
    """Constraint object for fixing a bond length."""
    def __init__(self, a1, a2):
        """Fix distance between atoms with indices a1 and a2."""
        self.indices = [a1, a2]

    def adjust_positions(self, old, new):
        p1, p2 = old[self.indices]
        d = p2 - p1
        p = sqrt(np.dot(d, d))
        q1, q2 = new[self.indices]
        d = q2 - q1
        q = sqrt(np.dot(d, d))
        d *= 0.5 * (p - q) / q
        new[self.indices] = (q1 - d, q2 + d)

    def adjust_forces(self, positions, forces):
        d = np.subtract.reduce(positions[self.indices])
        d2 = np.dot(d, d)
        d *= 0.5 * np.dot(np.subtract.reduce(forces[self.indices]), d) / d2
        forces[self.indices] += (-d, d)

    def index_shuffle(self, ind):
        'Shuffle the indices of the two atoms in this constraint'
        newa = [-1, -1] # Signal error
        for new, old in slice2enlist(ind):
            for i, a in enumerate(self.indices):
                if old == a:
                    newa[i] = new
        if newa[0] == -1 or newa[1] == -1:
            raise IndexError('Constraint not part of slice')
        self.indices = newa

    def copy(self):
        return FixBondLength(*self.indices)

    def __repr__(self):
        return 'FixBondLength(%d, %d)' % tuple(self.indices)

class FixedMode(FixConstraint):
    """Constrain atoms to move along directions orthogonal to
    a given mode only."""

    def __init__(self, mode):
        self.mode = (np.asarray(mode) / np.sqrt((mode **2).sum())).reshape(-1)

    def adjust_positions(self, oldpositions, newpositions):
        newpositions = newpositions.ravel()
        oldpositions = oldpositions.ravel()
        step = newpositions - oldpositions
        newpositions -= self.mode * np.dot(step, self.mode)
        newpositions = newpositions.reshape(-1, 3)
        oldpositions = oldpositions.reshape(-1, 3)

    def adjust_forces(self, positions, forces):
        forces = forces.ravel()
        forces -= self.mode * np.dot(forces, self.mode)
        forces = forces.reshape(-1, 3)

    def index_shuffle(self, ind):
        eps = 1e-12
        mode = self.mode.reshape(-1, 3)
        excluded = np.ones(len(mode), dtype=bool)
        excluded[ind] = False
        if (abs(mode[excluded]) > eps).any():
            raise IndexError('All nonzero parts of mode not in slice')
        self.mode = mode[ind].ravel()

    def copy(self):
        return FixedMode(self.mode)

    def __repr__(self):
        return 'FixedMode(%s)' % self.mode.tolist()

class FixedPlane(FixConstraintSingle):
    """Constrain an atom *a* to move in a given plane only.

    The plane is defined by its normal: *direction*."""

    def __init__(self, a, direction):
        self.a = a
        self.dir = np.asarray(direction) / sqrt(np.dot(direction, direction))

    def adjust_positions(self, oldpositions, newpositions):
        step = newpositions[self.a] - oldpositions[self.a]
        newpositions[self.a] -= self.dir * np.dot(step, self.dir)

    def adjust_forces(self, positions, forces):
        forces[self.a] -= self.dir * np.dot(forces[self.a], self.dir)

    def copy(self):
        return FixedPlane(self.a, self.dir)

    def __repr__(self):
        return 'FixedPlane(%d, %s)' % (self.a, self.dir.tolist())


class FixedLine(FixConstraintSingle):
    """Constrain an atom *a* to move on a given line only.

    The line is defined by its *direction*."""

    def __init__(self, a, direction):
        self.a = a
        self.dir = np.asarray(direction) / sqrt(np.dot(direction, direction))

    def adjust_positions(self, oldpositions, newpositions):
        step = newpositions[self.a] - oldpositions[self.a]
        x = np.dot(step, self.dir)
        newpositions[self.a] = oldpositions[self.a] + x * self.dir

    def adjust_forces(self, positions, forces):
        forces[self.a] = self.dir * np.dot(forces[self.a], self.dir)

    def copy(self):
        return FixedLine(self.a, self.dir)

    def __repr__(self):
        return 'FixedLine(%d, %s)' % (self.a, self.dir.tolist())

class FixCartesian(FixConstraintSingle):
    "Fix an atom in the directions of the cartesian coordinates."
    def __init__(self, a, mask=(1, 1, 1)):
        self.a = a
        self.mask = -(np.array(mask) - 1)

    def adjust_positions(self, old, new):
        step = new[self.a] - old[self.a]
        step *= self.mask
        new[self.a] = old[self.a] + step

    def adjust_forces(self, positions, forces):
        forces[self.a] *= self.mask

    def copy(self):
        return FixCartesian(self.a, 1 - self.mask)

    def __repr__(self):
        return 'FixCartesian(indice=%s mask=%s)' % (self.a, self.mask)

class fix_cartesian(FixCartesian):
    "Backwards compatibility for FixCartesian."
    def __init__(self, a, mask=(1, 1, 1)):
        import warnings
        super(fix_cartesian, self).__init__(a, mask)
        warnings.warn('fix_cartesian is deprecated. Please use FixCartesian'
                      ' instead.', DeprecationWarning, stacklevel=2)

class FixScaled(FixConstraintSingle):
    "Fix an atom in the directions of the unit vectors."
    def __init__(self, cell, a, mask=(1, 1, 1)):
        self.cell = cell
        self.a = a
        self.mask = np.array(mask)

    def adjust_positions(self, old, new):
        scaled_old = np.linalg.solve(self.cell.T, old.T).T
        scaled_new = np.linalg.solve(self.cell.T, new.T).T
        for n in range(3):
            if self.mask[n]:
                scaled_new[self.a, n] = scaled_old[self.a, n]
        new[self.a] = np.dot(scaled_new, self.cell)[self.a]

    def adjust_forces(self, positions, forces):
        scaled_forces = np.linalg.solve(self.cell.T, forces.T).T
        scaled_forces[self.a] *= -(self.mask - 1)
        forces[self.a] = np.dot(scaled_forces, self.cell)[self.a]

    def copy(self):
        return FixScaled(self.cell, self.a, self.mask)

    def __repr__(self):
        return 'FixScaled(%s, %d, %s)' % (repr(self.cell),
                                          self.a,
                                          repr(self.mask))

class fix_scaled(FixScaled):
    "Backwards compatibility for FixScaled."
    def __init__(self, cell, a, mask=(1, 1, 1)):
        import warnings
        super(fix_scaled, self).__init__(cell, a, mask)
        warnings.warn('fix_scaled is deprecated. Please use FixScaled '
                      'instead.', DeprecationWarning, stacklevel=2)

class FixInternals(FixConstraint):
    """Constraint object for fixing multiple internal coordinates.

    Allows fixing bonds, angles, and dihedrals."""
    def __init__(self, atoms=None, bonds=None, angles=None, dihedrals=None,
             epsilon=1.e-7, _copy_init=None):
        if _copy_init is None:
            if atoms is None:
                raise ValueError('Atoms object has to be defined.')
            masses = atoms.get_masses()
            if bonds is None:
                bonds = []
            if angles is None:
                angles = []
            if dihedrals is None:
                dihedrals = []
            self.n = len(bonds) + len(angles) + len(dihedrals)
            self.constraints = []
            for bond in bonds:
                masses_bond = masses[bond[1]]
                self.constraints.append(self.FixBondLengthAlt(bond[0], bond[1],
                                                              masses_bond))
            for angle in angles:
                masses_angle = masses[angle[1]]
                self.constraints.append(self.FixAngle(angle[0], angle[1],
                                                      masses_angle))
            for dihedral in dihedrals:
                masses_dihedral = masses[dihedral[1]]
                self.constraints.append(self.FixDihedral(dihedral[0],
                                                         dihedral[1],
                                                         masses_dihedral))
            self.epsilon = epsilon

        #copy case for __init__
        else:
            self.constraints = _copy_init
            self.n = len(self.constraints)
            self.epsilon = epsilon

    def adjust_positions(self, old, new):
        for constraint in self.constraints:
            constraint.set_h_vectors(old)
        for j in range(50):
            maxerr = 0.0
            for constraint in self.constraints:
                constraint.adjust_positions(old, new)
                maxerr = max(abs(constraint.sigma), maxerr)
            if maxerr < self.epsilon:
                return
        raise ValueError('Shake did not converge.')

    def adjust_forces(self, positions, forces):
        #Project out translations and rotations and all other constraints
        N = len(forces)
        list2_constraints = list(np.zeros((6, N, 3)))
        tx, ty, tz, rx, ry, rz = list2_constraints

        list_constraints = [r.ravel() for r in list2_constraints]
        
        tx[:, 0] = 1.0
        ty[:, 1] = 1.0
        tz[:, 2] = 1.0
        ff = forces.ravel()
        
        #Calculate the center of mass
        center = positions.sum(axis=0) / N

        rx[:, 1] = -(positions[:, 2] - center[2])
        rx[:, 2] = positions[:, 1] - center[1]
        ry[:, 0] = positions[:, 2] - center[2]
        ry[:, 2] = -(positions[:, 0] - center[0])
        rz[:, 0] = -(positions[:, 1] - center[1])
        rz[:, 1] = positions[:, 0] - center[0]
        
        #Normalizing transl., rotat. constraints
        for r in list2_constraints:
            r /= np.linalg.norm(r.ravel())
        
        #Add all angle, etc. constraint vectors
        for constraint in self.constraints:
            constraint.adjust_forces(positions, forces)
            list_constraints.insert(0, constraint.h)
        #QR DECOMPOSITION - GRAM SCHMIDT

        list_constraints = [r.ravel() for r in list_constraints]
        aa = np.column_stack(list_constraints)
        (aa, bb) = np.linalg.qr(aa, mode = 'full')
        #Projektion
        hh = []
        for i, constraint in enumerate(self.constraints):
            hh.append(aa[:, i] * np.row_stack(aa[:, i]))

        txx = aa[:, self.n] * np.row_stack(aa[:, self.n])
        tyy = aa[:, self.n+1] * np.row_stack(aa[:, self.n+1])
        tzz = aa[:, self.n+2] * np.row_stack(aa[:, self.n+2])
        rxx = aa[:, self.n+3] * np.row_stack(aa[:, self.n+3])
        ryy = aa[:, self.n+4] * np.row_stack(aa[:, self.n+4])
        rzz = aa[:, self.n+5] * np.row_stack(aa[:, self.n+5])
        T = txx + tyy + tzz + rxx + ryy + rzz
        for vec in hh:
            T += vec
        ff = np.dot(T, np.row_stack(ff))
        forces[:, :] -= np.dot(T, np.row_stack(ff)).reshape(-1, 3)

    def copy(self):
        return FixInternals(epsilon=self.epsilon, _copy_init=self.constraints)

    def __repr__(self):
        constraints = repr(self.constraints)
        return 'FixInternals(_copy_init=%s, epsilon=%s)' % (constraints,
                                                            repr(self.epsilon))
    
    def __str__(self):
        return '\n'.join([repr(c) for c in self.constraints])

    #Classes for internal use in FixInternals
    class FixBondLengthAlt:
        """Constraint subobject for fixing bond length within FixInternals."""
        def __init__(self, bond, indices, masses, maxstep=0.01):
            """Fix distance between atoms with indices a1, a2."""
            self.indices = indices
            self.bond = bond
            self.h1, self.h2 = None
            self.masses = masses
            self.h = []
            self.sigma = 1.

        def set_h_vectors(self, pos):
            dist1 = pos[self.indices[0]] - pos[self.indices[1]]
            self.h1 = 2 * dist1
            self.h2 = -self.h1

        def adjust_positions(self, old, new):
            h1 = self.h1 / self.masses[0]
            h2 = self.h2 / self.masses[1]
            dist1 = new[self.indices[0]] - new[self.indices[1]]
            dist = np.dot(dist1, dist1)
            self.sigma = dist - self.bond**2
            lamda = -self.sigma / (2 * np.dot(dist1, (h1 - h2)))
            new[self.indices[0]] += lamda * h1
            new[self.indices[1]] += lamda * h2

        def adjust_forces(self, positions, forces):
            self.h1 = 2 * (positions[self.indices[0]] -
                           positions[self.indices[1]])
            self.h2 = -self.h1
            self.h = np.zeros([len(forces)*3])
            self.h[(self.indices[0])*3]   = self.h1[0]
            self.h[(self.indices[0])*3+1] = self.h1[1]
            self.h[(self.indices[0])*3+2] = self.h1[2]
            self.h[(self.indices[1])*3]   = self.h2[0]
            self.h[(self.indices[1])*3+1] = self.h2[1]
            self.h[(self.indices[1])*3+2] = self.h2[2]
            self.h /= np.linalg.norm(self.h)

        def __repr__(self):
            return 'FixBondLengthAlt(%d, %d, %d)' % \
                tuple(self.bond, self.indices)

    class FixAngle:
        """Constraint object for fixing an angle within
        FixInternals."""
        def __init__(self, angle, indices, masses):
            """Fix atom movement to construct a constant angle."""
            self.indices = indices
            self.a1m, self.a2m, self.a3m = masses
            self.angle = np.cos(angle)
            self.h1 = self.h2 = self.h3 = None
            self.h = []
            self.sigma = 1.

        def set_h_vectors(self, pos):
            r21 = pos[self.indices[0]] - pos[self.indices[1]]
            r21_len = np.linalg.norm(r21)
            e21 = r21 / r21_len
            r23 = pos[self.indices[2]] - pos[self.indices[1]]
            r23_len = np.linalg.norm(r23)
            e23 = r23 / r23_len
            angle = np.dot(e21, e23)
            self.h1 = -2 * angle * ((angle * e21 - e23) / (r21_len))
            self.h3 = -2 * angle * ((angle * e23 - e21) / (r23_len))
            self.h2 = -(self.h1 + self.h3)

        def adjust_positions(self, oldpositions, newpositions):
            r21 = newpositions[self.indices[0]] - newpositions[self.indices[1]]
            r21_len = np.linalg.norm(r21)
            e21 = r21 / r21_len
            r23 = newpositions[self.indices[2]] - newpositions[self.indices[1]]
            r23_len = np.linalg.norm(r23)
            e23 = r23 / r23_len
            angle = np.dot(e21, e23)
            self.sigma = (angle - self.angle) * (angle + self.angle)
            h1 = self.h1 / self.a1m
            h3 = self.h3 / self.a3m
            h2 = self.h2 / self.a2m
            h21 = h1 - h2
            h23 = h3 - h2
            # Calculating new positions
            deriv = (((np.dot(r21, h23) + np.dot(r23, h21))
                      / (r21_len * r23_len))
                     - (np.dot(r21, h21) / (r21_len * r21_len)
                        + np.dot(r23, h23) / (r23_len * r23_len)) * angle)
            deriv *= 2 * angle
            lamda = -self.sigma / deriv
            newpositions[self.indices[0]] += lamda * h1
            newpositions[self.indices[1]] += lamda * h2
            newpositions[self.indices[2]] += lamda * h3

        def adjust_forces(self, positions, forces):
            r21 = positions[self.indices[0]] - positions[self.indices[1]]
            r21_len = np.linalg.norm(r21)
            e21 = r21 / r21_len
            r23 = positions[self.indices[2]] - positions[self.indices[1]]
            r23_len = np.linalg.norm(r23)
            e23 = r23 / r23_len
            angle = np.dot(e21, e23)
            self.h1 = -2 * angle * (angle * e21 - e23) / r21_len
            self.h3 = -2 * angle * (angle * e23 - e21) / r23_len
            self.h2 = -(self.h1 + self.h3)
            self.h = np.zeros([len(positions)*3])
            self.h[(self.indices[0])*3]   = self.h1[0]
            self.h[(self.indices[0])*3+1] = self.h1[1]
            self.h[(self.indices[0])*3+2] = self.h1[2]
            self.h[(self.indices[1])*3]   = self.h2[0]
            self.h[(self.indices[1])*3+1] = self.h2[1]
            self.h[(self.indices[1])*3+2] = self.h2[2]
            self.h[(self.indices[2])*3]   = self.h3[0]
            self.h[(self.indices[2])*3+1] = self.h3[1]
            self.h[(self.indices[2])*3+2] = self.h3[2]
            self.h /= np.linalg.norm(self.h)

        def __repr__(self):
            return 'FixAngle(%s, %f)' % (tuple(self.indices),
                                         np.arccos(self.angle))
                

    class FixDihedral:
        """Constraint object for fixing an dihedral using
        the shake algorithm. This one allows also other constraints."""
        def __init__(self, angle, indices, masses):
            """Fix atom movement to construct a constant dihedral angle."""
            self.indices = indices
            self.a1m, self.a2m, self.a3m, self.a4m = masses
            self.angle = np.cos(angle)
            self.h1 = self.h2 = self.h3 = self.h4 = None
            self.h = []
            self.sigma = 1.

        def set_h_vectors(self, pos):
            r12 = pos[self.indices[1]] - pos[self.indices[0]]
            r12_len = np.linalg.norm(r12)
            e12 = r12 / r12_len
            r23 = pos[self.indices[2]] - pos[self.indices[1]]
            r23_len = np.linalg.norm(r23)
            e23 = r23 / r23_len
            r34 = pos[self.indices[3]] - pos[self.indices[2]]
            r34_len = np.linalg.norm(r34)
            e34 = r34 / r34_len
            a = -r12 - np.dot(-r12, e23) * e23
            a_len = np.linalg.norm(a)
            ea = a / a_len
            b = r34 - np.dot(r34, e23) * e23
            b_len = np.linalg.norm(b)
            eb = b / b_len
            angle = np.dot(ea, eb).clip(-1.0, 1.0)
            self.h1 = (eb - angle * ea) / a_len
            self.h4 = (ea - angle * eb) / b_len
            self.h2 = self.h1 * (np.dot(-r12, e23) / r23_len -1)
            self.h2 += np.dot(r34, e23) / r23_len * self.h4
            self.h3 = -self.h4 * (np.dot(r34, e23) / r23_len + 1)
            self.h3 += np.dot(r12, e23) / r23_len * self.h1

        def adjust_positions(self, oldpositions, newpositions):
            r12 = newpositions[self.indices[1]] - newpositions[self.indices[0]]
            r12_len = np.linalg.norm(r12)
            e12 = r12 / r12_len
            r23 = newpositions[self.indices[2]] - newpositions[self.indices[1]]
            r23_len = np.linalg.norm(r23)
            e23 = r23 / r23_len
            r34 = newpositions[self.indices[3]] - newpositions[self.indices[2]]
            r34_len = np.linalg.norm(r34)
            e34 = r34 / r34_len
            n1 = np.cross(r12, r23)
            n1_len = np.linalg.norm(n1)
            n1e = n1 / n1_len
            n2 = np.cross(r23, r34)
            n2_len = np.linalg.norm(n2)
            n2e = n2 / n2_len
            angle = np.dot(n1e, n2e).clip(-1.0, 1.0)
            self.sigma = (angle - self.angle) * (angle + self.angle)
            h1 = self.h1 / self.a1m
            h2 = self.h2 / self.a2m
            h3 = self.h3 / self.a3m
            h4 = self.h4 / self.a4m
            h12 = h2 - h1
            h23 = h3 - h2
            h34 = h4 - h3
            deriv = ((np.dot(n1, np.cross(r34, h23) + np.cross(h34, r23))
                      + np.dot(n2, np.cross(r23, h12) + np.cross(h23, r12)))
                     / (n1_len * n2_len))
            deriv -= (((np.dot(n1, np.cross(r23, h12) + np.cross(h23, r12))
                        / n1_len**2)
                       + (np.dot(n2, np.cross(r34, h23) + np.cross(h34, r23))
                          / n2_len**2)) * angle)
            deriv *= -2 * angle
            lamda = -self.sigma / deriv
            newpositions[self.indices[0]] += lamda * h1
            newpositions[self.indices[1]] += lamda * h2
            newpositions[self.indices[2]] += lamda * h3
            newpositions[self.indices[3]] += lamda * h4

        def adjust_forces(self, positions, forces):
            r12 = positions[self.indices[1]] - positions[self.indices[0]]
            r12_len = np.linalg.norm(r12)
            e12 = r12 / r12_len
            r23 = positions[self.indices[2]] - positions[self.indices[1]]
            r23_len = np.linalg.norm(r23)
            e23 = r23 / r23_len
            r34 = positions[self.indices[3]] - positions[self.indices[2]]
            r34_len = np.linalg.norm(r34)
            e34 = r34 / r34_len
            a = -r12 - np.dot(-r12, e23) * e23
            a_len = np.linalg.norm(a)
            ea = a / a_len
            b = r34 - np.dot(r34, e23) * e23
            b_len = np.linalg.norm(b)
            eb = b / b_len
            angle = np.dot(ea, eb).clip(-1.0, 1.0)
            self.h1 = (eb - angle * ea) / a_len
            self.h4 = (ea - angle * eb) / b_len
            self.h2 = self.h1 * (np.dot(-r12, e23) / r23_len - 1)
            self.h2 += np.dot(r34, e23) / r23_len * self.h4
            self.h3 = -self.h4 * (np.dot(r34, e23) / r23_len + 1)
            self.h3 -= np.dot(-r12, e23) / r23_len * self.h1

            self.h = np.zeros([len(positions)*3])
            self.h[(self.indices[0])*3]   = self.h1[0]
            self.h[(self.indices[0])*3+1] = self.h1[1]
            self.h[(self.indices[0])*3+2] = self.h1[2]
            self.h[(self.indices[1])*3]   = self.h2[0]
            self.h[(self.indices[1])*3+1] = self.h2[1]
            self.h[(self.indices[1])*3+2] = self.h2[2]
            self.h[(self.indices[2])*3]   = self.h3[0]
            self.h[(self.indices[2])*3+1] = self.h3[1]
            self.h[(self.indices[2])*3+2] = self.h3[2]
            self.h[(self.indices[3])*3]   = self.h4[0]
            self.h[(self.indices[3])*3+1] = self.h4[1]
            self.h[(self.indices[3])*3+2] = self.h4[2]
            self.h /= np.linalg.norm(self.h)
        
        def __repr__(self):
            return 'FixDihedral(%s, %f)' % (tuple(self.indices), self.angle)

class BondSpring(FixConstraint):
    """Forces two atoms to stay close together by applying no force if they
    are below threshhold_length, and applying a Hookian force when the
    distance between them exceeds the thresshhold_length.
    
    a1, a2 : indices of atoms 1 and 2
    a2 can alternately be a position in space to tether a1 to
    threshhold_length (float) : the length below which there is no force
    springconstant (integer) : Hook's law constant to apply when distance
        between the two atoms exceeds threshhold_length, dimensions of 
        (force / length)
    """
    def __init__(self, a1, a2, threshhold_length, springconstant):
        if type(a2) == int:
            self._type = 2 # two atoms tethered together
            self.indices = [a1, a2]
        else:
            self._type = 1 # one atom tethered to a point in space
            self.index = a1
            self.origin = np.array(a2)
        self.threshhold = threshhold_length
        self.spring = springconstant

    def adjust_positions(self, oldpositions, newpositions):
        pass

    def adjust_forces(self, positions, forces):
        if self._type == 2:
            p1, p2 = positions[self.indices]
        else:
            p1 = positions[self.index]
            p2 = self.origin
        displace = p2 - p1
        bondlength = np.linalg.norm(displace)
        if bondlength > self.threshhold:
            magnitude = self.spring * (bondlength - self.threshhold)
            direction = displace / np.linalg.norm(displace)
            if self._type == 2:
                forces[self.indices[0]] += direction * magnitude
                forces[self.indices[1]] -= direction * magnitude
            else:
                forces[self.index] += direction * magnitude

    def __repr__(self):
        if self._type == 2:
            return 'BondSpring(%d, %d)' % tuple(self.indices)
        else:
            return 'BondSpring(%d) to cartesian' % self.index

    def copy(self):
        if self._type == 2:
            return BondSpring(a1=self.indices[0], a2=self.indices[1],
                              threshhold_length=self.threshhold,
                              springconstant=self.spring)
        else:
            return BondSpring(a1=self.index, a2=self.origin,
                              threshhold_length=self.threshhold,
                              springconstant=self.spring)


class Filter:
    """Subset filter class."""
    def __init__(self, atoms, indices=None, mask=None):
        """Filter atoms.

        This filter can be used to hide degrees of freedom in an Atoms
        object.

        Parameters
        ----------
        indices : list of int
           Indices for those atoms that should remain visible.
        mask : list of bool
           One boolean per atom indicating if the atom should remain
           visible or not.
        """

        self.atoms = atoms
        self.constraints = []

        if indices is None and mask is None:
            raise ValueError('Use "indices" or "mask".')
        if indices is not None and mask is not None:
            raise ValueError('Use only one of "indices" and "mask".')

        if mask is not None:
            self.index = np.asarray(mask, bool)
            self.n = self.index.sum()
        else:
            self.index = np.asarray(indices, int)
            self.n = len(self.index)

    def get_cell(self):
        """Returns the computational cell.

        The computational cell is the same as for the original system.
        """
        return self.atoms.get_cell()

    def get_pbc(self):
        """Returns the periodic boundary conditions.

        The boundary conditions are the same as for the original system.
        """
        return self.atoms.get_pbc()

    def get_positions(self):
        "Return the positions of the visible atoms."
        return self.atoms.get_positions()[self.index]

    def set_positions(self, positions):
        "Set the positions of the visible atoms."
        pos = self.atoms.get_positions()
        pos[self.index] = positions
        self.atoms.set_positions(pos)

    positions = property(get_positions, set_positions, 
                         doc='Positions of the atoms')

    def get_momenta(self):
        "Return the momenta of the visible atoms."
        return self.atoms.get_momenta()[self.index]

    def set_momenta(self, momenta):
        "Set the momenta of the visible atoms."
        mom = self.atoms.get_momenta()
        mom[self.index] = momenta
        self.atoms.set_momenta(mom)

    def get_atomic_numbers(self):
        "Return the atomic numbers of the visible atoms."
        return self.atoms.get_atomic_numbers()[self.index]

    def set_atomic_numbers(self, atomic_numbers):
        "Set the atomic numbers of the visible atoms."
        z = self.atoms.get_atomic_numbers()
        z[self.index] = atomic_numbers
        self.atoms.set_atomic_numbers(z)

    def get_tags(self):
        "Return the tags of the visible atoms."
        return self.atoms.get_tags()[self.index]

    def set_tags(self, tags):
        "Set the tags of the visible atoms."
        tg = self.atoms.get_tags()
        tg[self.index] = tags
        self.atoms.set_tags(tg)

    def get_forces(self, *args, **kwargs):
        return self.atoms.get_forces(*args, **kwargs)[self.index]

    def get_stress(self):
        return self.atoms.get_stress()

    def get_stresses(self):
        return self.atoms.get_stresses()[self.index]

    def get_masses(self):
        return self.atoms.get_masses()[self.index]

    def get_potential_energy(self):
        """Calculate potential energy.

        Returns the potential energy of the full system.
        """
        return self.atoms.get_potential_energy()

    def get_chemical_symbols(self):
        return self.atoms.get_chemical_symbols()

    def get_initial_magnetic_moments(self):
        return self.atoms.get_initial_magnetic_moments()

    def get_calculator(self):
        """Returns the calculator.

        WARNING: The calculator is unaware of this filter, and sees a
        different number of atoms.
        """
        return self.atoms.get_calculator()

    def has(self, name):
        """Check for existance of array."""
        return self.atoms.has(name)

    def __len__(self):
        "Return the number of movable atoms."
        return self.n

    def __getitem__(self, i):
        "Return an atom."
        return self.atoms[self.index[i]]


class StrainFilter(Filter):
    """Modify the supercell while keeping the scaled positions fixed.

    Presents the strain of the supercell as the generalized positions,
    and the global stress tensor (times the volume) as the generalized
    force.

    This filter can be used to relax the unit cell until the stress is
    zero.  If MDMin is used for this, the timestep (dt) to be used
    depends on the system size. 0.01/x where x is a typical dimension
    seems like a good choice.

    The stress and strain are presented as 6-vectors, the order of the
    components follow the standard engingeering practice: xx, yy, zz,
    yz, xz, xy.

    """
    def __init__(self, atoms, mask=None):
        """Create a filter applying a homogeneous strain to a list of atoms.

        The first argument, atoms, is the atoms object.

        The optional second argument, mask, is a list of six booleans,
        indicating which of the six independent components of the
        strain that are allowed to become non-zero.  It defaults to
        [1,1,1,1,1,1].

        """

        self.atoms = atoms
        self.strain = np.zeros(6)

        if mask is None:
            self.mask = np.ones(6)
        else:
            self.mask = np.array(mask)

        self.index = np.arange(len(atoms))
        self.n = self.index.sum()

        self.origcell = atoms.get_cell()

    def get_positions(self):
        return self.strain.reshape((2, 3))

    def set_positions(self, new):
        new = new.ravel() * self.mask
        eps = np.array([[1.0 + new[0], 0.5 * new[5], 0.5 * new[4]],
                        [0.5 * new[5], 1.0 + new[1], 0.5 * new[3]],
                        [0.5 * new[4], 0.5 * new[3], 1.0 + new[2]]])

        self.atoms.set_cell(np.dot(self.origcell, eps), scale_atoms=True)
        self.strain[:] = new

    def get_forces(self):
        stress = self.atoms.get_stress()
        return -self.atoms.get_volume() * (stress * self.mask).reshape((2, 3))

    def get_potential_energy(self):
        return self.atoms.get_potential_energy()

    def has(self, x):
        return self.atoms.has(x)

    def __len__(self):
        return 2

class UnitCellFilter(Filter):
    """Modify the supercell and the atom positions. """
    def __init__(self, atoms, mask=None):
        """Create a filter that returns the atomic forces and unit
        cell stresses together, so they can simultaneously be
        minimized.

        The first argument, atoms, is the atoms object.

        The optional second argument, mask, is a list of booleans,
        indicating which of the six independent
        components of the strain are relaxed.
        1, True = relax to zero
        0, False = fixed, ignore this component

        use atom Constraints, e.g. FixAtoms, to control relaxation of
        the atoms.

        #this should be equivalent to the StrainFilter
        >>> atoms = Atoms(...)
        >>> atoms.set_constraint(FixAtoms(mask=[True for atom in atoms]))
        >>> ucf = UCFilter(atoms)

        You should not attach this UCFilter object to a
        trajectory. Instead, create a trajectory for the atoms, and
        attach it to an optimizer like this:

        >>> atoms = Atoms(...)
        >>> ucf = UCFilter(atoms)
        >>> qn = QuasiNewton(ucf)
        >>> traj = PickleTrajectory('TiO2.traj','w',atoms)
        >>> qn.attach(traj)
        >>> qn.run(fmax=0.05)

        Helpful conversion table
        ========================
        0.05 eV/A^3   = 8 GPA
        0.003 eV/A^3  = 0.48 GPa
        0.0006 eV/A^3 = 0.096 GPa
        0.0003 eV/A^3 = 0.048 GPa
        0.0001 eV/A^3 = 0.02 GPa
        """

        Filter.__init__(self, atoms, indices=range(len(atoms)))

        self.atoms = atoms
        self.strain = np.zeros(6)

        if mask is None:
            self.mask = np.ones(6)
        else:
            self.mask = np.array(mask)

        self.origcell = atoms.get_cell()

    def get_positions(self):
        '''
        this returns an array with shape (natoms + 2,3).

        the first natoms rows are the positions of the atoms, the last
        two rows are the strains associated with the unit cell
        '''

        atom_positions = self.atoms.get_positions()
        strains = self.strain.reshape((2, 3))

        natoms = len(self.atoms)
        all_pos = np.zeros((natoms + 2, 3), np.float)
        all_pos[0:natoms, :] = atom_positions
        all_pos[natoms:, :] = strains

        return all_pos

    def set_positions(self, new):
        '''
        new is an array with shape (natoms+2,3).

        the first natoms rows are the positions of the atoms, the last
        two rows are the strains used to change the cell shape.

        The atom positions are set first, then the unit cell is
        changed keeping the atoms in their scaled positions.
        '''

        natoms = len(self.atoms)

        atom_positions = new[0:natoms, :]
        self.atoms.set_positions(atom_positions)

        new = new[natoms:, :] #this is only the strains
        new = new.ravel() * self.mask
        eps = np.array([[1.0 + new[0], 0.5 * new[5], 0.5 * new[4]],
                        [0.5 * new[5], 1.0 + new[1], 0.5 * new[3]],
                        [0.5 * new[4], 0.5 * new[3], 1.0 + new[2]]])

        self.atoms.set_cell(np.dot(self.origcell, eps), scale_atoms=True)
        self.strain[:] = new

    def get_forces(self, apply_constraint=False):
        '''
        returns an array with shape (natoms+2,3) of the atomic forces
        and unit cell stresses.

        the first natoms rows are the forces on the atoms, the last
        two rows are the stresses on the unit cell, which have been
        reshaped to look like "atomic forces". i.e.,

        f[-2] = -vol*[sxx,syy,szz]*mask[0:3]
        f[-1] = -vol*[syz, sxz, sxy]*mask[3:]

        apply_constraint is an argument expected by ase
        '''

        stress = self.atoms.get_stress()
        atom_forces = self.atoms.get_forces()

        natoms = len(self.atoms)
        all_forces = np.zeros((natoms+2, 3), np.float)
        all_forces[0:natoms, :] = atom_forces

        vol = self.atoms.get_volume()
        stress_forces = -vol * (stress * self.mask).reshape((2, 3))
        all_forces[natoms:, :] = stress_forces
        return all_forces

    def get_potential_energy(self):
        return self.atoms.get_potential_energy()

    def has(self, x):
        return self.atoms.has(x)

    def __len__(self):
        return (2 + len(self.atoms))
