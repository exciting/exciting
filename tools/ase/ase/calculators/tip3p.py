"""TIP3P potential, constraints and dynamics."""
from math import pi, sin, cos

import numpy as np

import ase.units as units
from ase.parallel import world
from ase.md.md import MolecularDynamics

qH = 0.417
sigma0 = 3.15061
epsilon0 = 0.1521 * units.kcal / units.mol
rOH = 0.9572
thetaHOH = 104.52 / 180 * pi


class TIP3P:
    def __init__(self, rc=9.0, width=1.0):
        self.energy = None
        self.forces = None

        self.rc1 = rc - width
        self.rc2 = rc
        
    def get_spin_polarized(self):
        return False
    
    def update(self, atoms):
        if (self.energy is None or
            len(self.numbers) != len(atoms) or
            (self.numbers != atoms.get_atomic_numbers()).any()):
            self.calculate(atoms)
        elif ((self.positions != atoms.get_positions()).any() or
              (self.pbc != atoms.get_pbc()).any() or
              (self.cell != atoms.get_cell()).any()):
            self.calculate(atoms)

    def calculation_required(self, atoms, quantities):
        if len(quantities) == 0:
            return False

        return (self.energy is None or
                len(self.numbers) != len(atoms) or
                (self.numbers != atoms.get_atomic_numbers()).any() or
                (self.positions != atoms.get_positions()).any() or
                (self.pbc != atoms.get_pbc()).any() or
                (self.cell != atoms.get_cell()).any())
                
    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()
    
    def get_stress(self, atoms):
        raise NotImplementedError
    
    def calculate(self, atoms):
        self.positions = atoms.get_positions().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()
        natoms = len(atoms)
        nH2O = natoms // 3

        assert self.pbc.all()
        C = self.cell.diagonal()
        assert not (self.cell - np.diag(C)).any()
        assert (C >= 2 * self.rc2).all()
        self.numbers = atoms.get_atomic_numbers()
        Z = self.numbers.reshape((-1, 3))
        assert (Z[:, 1:] == 1).all() and (Z[:, 0] == 8).all()

        R = self.positions.reshape((nH2O, 3, 3))
        RO = R[:, 0]
        
        self.energy = 0.0
        self.forces = np.zeros((natoms, 3))
        
        if world.size == 1:
            mya = range(nH2O - 1)
        else:
            rank = world.rank
            size = world.size
            assert nH2O // (2 * size) == 0
            mynH2O = nH2O // 2 // size
            mya = (range(rank * n, (rank + 1) * n) +
                   range((size - rank - 1) * n, (size - rank) * n))

        q = np.empty(3)
        q[:] = qH * (units.Hartree * units.Bohr)**0.5
        q[0] *= -2
        
        for a in mya:
            DOO = (RO[a + 1:] - RO[a] + 0.5 * C) % C - 0.5 * C
            dOO = (DOO**2).sum(axis=1)**0.5
            x1 = dOO > self.rc1
            x2 = dOO < self.rc2
            f = np.zeros(nH2O - a - 1)
            f[x2] = 1.0
            dfdd = np.zeros(nH2O - a - 1)
            x12 = np.logical_and(x1, x2)
            d = (dOO[x12] - self.rc1) / (self.rc2 - self.rc1)
            f[x12] -= d**2 * (3.0 - 2.0 * d)
            dfdd[x12] -= 6.0 / (self.rc2 - self.rc1) * d * (1.0 - d)

            y = (sigma0 / dOO)**6
            y2 = y**2
            e = 4 * epsilon0 * (y2 - y)
            self.energy += np.dot(e, f)
            dedd = 24 * epsilon0 * (2 * y2 - y) / dOO * f - e * dfdd
            F = (dedd / dOO)[:, np.newaxis] * DOO
            self.forces[(a + 1) * 3::3] += F
            self.forces[a * 3] -= F.sum(axis=0)

            for i in range(3):
                D = (R[a + 1:] - R[a, i] + 0.5 * C) % C - 0.5 * C
                d = (D**2).sum(axis=2)**0.5
                e = q[i] * q / d
                self.energy += np.dot(f, e).sum()
                F = (e / d**2 * f[:, np.newaxis])[:, :, np.newaxis] * D
                F[:, 0] -= (e.sum(axis=1) * dfdd / dOO)[:, np.newaxis] * DOO 
                self.forces[(a + 1) * 3:] += F.reshape((-1, 3))
                self.forces[a * 3 + i] -= F.sum(axis=0).sum(axis=0)

        self.energy = world.sum(self.energy)
        world.sum(self.forces)


class H2OConstraint:
    """Constraint object for a rigid H2O molecule."""
    def __init__(self, r=rOH, theta=thetaHOH, iterations=23, masses=None):
        self.r = r
        self.theta = theta
        self.iterations = iterations
        self.m = masses

    def set_masses(self, masses):
        self.m = masses

    def adjust_positions(self, old, new):
        bonds = [(0, 1, self.r), (0, 2, self.r)]
        if self.theta:
            bonds.append((1, 2, sin(self.theta / 2) * self.r * 2))
        for iter in range(self.iterations):
            for i, j, r in bonds:
                D = old[i::3] - old[j::3]
                m1 = self.m[i]
                m2 = self.m[j]
                a = new[i::3]
                b = new[j::3]
                B = a - b
                x = (D**2).sum(axis=1)
                y = (D * B).sum(axis=1)
                z = (B**2).sum(axis=1) - r**2
                k = m1 * m2 / (m1 + m2) * ((y**2 - x * z)**0.5 - y) / x
                k.shape = (-1, 1)
                a += k / m1 * D
                b -= k / m2 * D

    def adjust_forces(self, positions, forces):
        pass
        
    def copy(self):
        return H2OConstraint(self.r, self.theta, self.iterations, self.m)


class Verlet(MolecularDynamics):
    def step(self, f):
        atoms = self.atoms
        m = atoms.get_masses()[:, np.newaxis]
        v = self.atoms.get_velocities()
        r0 = atoms.get_positions()
        r = r0 + self.dt * v + self.dt**2 * f / m
        atoms.set_positions(r)
        r = atoms.get_positions()
        v = (r - r0) / self.dt
        self.atoms.set_velocities(v)
        return atoms.get_forces()
