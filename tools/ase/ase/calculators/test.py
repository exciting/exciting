from math import pi, ceil

import numpy as np

from ase.atoms import Atoms
from ase.parallel import world, rank, distribute_cpus
try:
    from gpaw.mpi import SerialCommunicator
except:
    pass

def make_test_dft_calculation():
    a = b = 2.0
    c = 6.0
    atoms = Atoms(positions=[(0, 0, c / 2)],
                  symbols='H',
                  pbc=(1, 1, 0),
                  cell=(a, b, c),
                  calculator=TestCalculator())
    return atoms


class TestCalculator:
    def __init__(self, nk=8):
        assert nk % 2 == 0
        bzk = []
        weights = []
        ibzk = []
        w = 1.0 / nk**2
        for i in range(-nk + 1, nk, 2):
            for j in range(-nk + 1, nk, 2):
                k = (0.5 * i / nk, 0.5 * j / nk, 0)
                bzk.append(k)
                if i >= j > 0:
                    ibzk.append(k)
                    if i == j:
                        weights.append(4 * w)
                    else:
                        weights.append(8 * w)
        assert abs(sum(weights) - 1.0) < 1e-12
        self.bzk = np.array(bzk)
        self.ibzk = np.array(ibzk)
        self.weights = np.array(weights)

        # Calculate eigenvalues and wave functions:
        self.init()

    def init(self):
        nibzk = len(self.weights)
        nbands = 1

        V = -1.0
        self.eps = 2 * V * (np.cos(2 * pi * self.ibzk[:, 0]) +
                            np.cos(2 * pi * self.ibzk[:, 1]))
        self.eps.shape = (nibzk, nbands)

        self.psi = np.zeros((nibzk, 20, 20, 60), complex)
        phi = np.empty((2, 2, 20, 20, 60))
        z = np.linspace(-1.5, 1.5, 60, endpoint=False)
        for i in range(2):
            x = np.linspace(0, 1, 20, endpoint=False) - i
            for j in range(2):
                y = np.linspace(0, 1, 20, endpoint=False) - j
                r = (((x[:, None]**2 +
                       y**2)[:, :, None] +
                      z**2)**0.5).clip(0, 1)
                phi = 1.0 - r**2 * (3.0 - 2.0 * r)
                phase = np.exp(pi * 2j * np.dot(self.ibzk, (i, j, 0)))
                self.psi += phase[:, None, None, None] * phi

    def get_pseudo_wave_function(self, band=0, kpt=0, spin=0):
        assert spin == 0 and band == 0
        return self.psi[kpt]

    def get_eigenvalues(self, kpt=0, spin=0):
        assert spin == 0
        return self.eps[kpt]

    def get_number_of_bands(self):
        return 1

    def get_k_point_weights(self):
        return self.weights

    def get_number_of_spins(self):
        return 1

    def get_fermi_level(self):
        return 0.0


class TestPotential:
    def get_forces(self, atoms):
        E = 0.0
        R = atoms.positions
        F = np.zeros_like(R)
        for a, r in enumerate(R):
            D = R - r
            d = (D**2).sum(1)**0.5
            x = d - 1.0
            E += np.vdot(x, x)
            d[a] = 1
            F -= (x / d)[:, None] * D
        self.energy = 0.25 * E
        return F

    def get_potential_energy(self, atoms):
        self.get_forces(atoms)
        return self.energy

    def get_stress(self, atoms):
        raise NotImplementedError


def numeric_force(atoms, a, i, d=0.001):
    """Evaluate force along i'th axis on a'th atom using finite difference.

    This will trigger two calls to get_potential_energy(), with atom a moved
    plus/minus d in the i'th axial direction, respectively.
    """
    p0 = atoms.positions[a, i]
    atoms.positions[a, i] += d
    eplus = atoms.get_potential_energy()
    atoms.positions[a, i] -= 2 * d
    eminus = atoms.get_potential_energy()
    atoms.positions[a, i] = p0
    return (eminus - eplus) / (2 * d)


def numeric_forces(atoms, indices=None, axes=(0, 1, 2), d=0.001,
                   parallel=None):
    """Evaluate finite-difference forces on several atoms.

    Returns an array of forces for each specified atomic index and
    each specified axis, calculated using finite difference on each
    atom and direction separately.  Array has same shape as if
    returned from atoms.get_forces(); uncalculated elements are zero.

    Calculates all forces by default."""

    if indices is None:
        indices = range(len(atoms))
    F_ai = np.zeros_like(atoms.positions)
    n = len(indices) * len(axes)
    if parallel is None:
        atom_tasks = [atoms] * n
        master = True
    else:
        calc_comm, tasks_comm, tasks_rank = distribute_cpus(parallel, world)
        master = calc_comm.rank == 0 
        calculator = atoms.get_calculator()
        calculator.set(communicator=calc_comm)
        atom_tasks = [None] * n
        for i in range(n):
            if ((i - tasks_rank) % tasks_comm.size) == 0:
                atom_tasks[i] = atoms
    for ia, a in enumerate(indices):
        for ii, i in enumerate(axes):
            atoms = atom_tasks[ia * len(axes) + ii]
            if atoms is not None:
                print '# rank', rank, 'calculating atom', a, 'xyz'[i]
                force = numeric_force(atoms, a, i, d)
                if master:
                    F_ai[a, i] = force
    if parallel is not None:
        world.sum(F_ai)
    return F_ai
