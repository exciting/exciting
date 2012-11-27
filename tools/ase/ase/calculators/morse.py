from math import exp, sqrt
import numpy as np

class MorsePotential:
    """Morse potential.

    Default values chosen to be similar as Lennard-Jones.
    """
    def __init__(self, rho0=6.0, epsilon=1.0, r0=1.0):
        self.epsilon = epsilon
        self.rho0 = rho0
        self.r0 = r0
        self.positions = None

    def update(self, atoms):
        assert not atoms.get_pbc().any()
        if (self.positions is None or
            (self.positions != atoms.get_positions()).any()):
            self.calculate(atoms)

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        return self._forces

    def get_stress(self, atoms):
        return np.zeros((3, 3))
    
    def calculate(self, atoms):
        positions = atoms.get_positions()
        self.energy = 0.0
        self._forces = np.zeros((len(atoms), 3))
        preF = 2 * self.epsilon * self.rho0 / self.r0
        for i1, p1 in enumerate(positions):
            for i2, p2 in enumerate(positions[:i1]):
                diff = p2 - p1
                r = sqrt(np.dot(diff, diff))
                expf = exp(self.rho0 * (1.0 - r / self.r0))
                self.energy += self.epsilon * expf * (expf - 2)
                F = preF * expf * (expf - 1) * diff / r
                self._forces[i1] -= F
                self._forces[i2] += F
        self.positions = positions.copy()
