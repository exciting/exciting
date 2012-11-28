import numpy as np


class LennardJones:
    def __init__(self, epsilon=1.0, sigma=1.0):
        self.epsilon = epsilon
        self.sigma = sigma
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
        for i1, p1 in enumerate(positions):
            for i2, p2 in enumerate(positions[:i1]):
                diff = p2 - p1
                d2 = np.dot(diff, diff)
                c6 = (self.sigma**2 / d2)**3
                c12 = c6**2
                self.energy += 4 * self.epsilon * (c12 - c6)
                F = 24 * self.epsilon * (2 * c12 - c6) / d2 * diff
                self._forces[i1] -= F
                self._forces[i2] += F
        self.positions = positions.copy()
