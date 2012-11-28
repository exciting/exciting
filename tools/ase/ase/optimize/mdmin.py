import numpy as np

from ase.optimize.optimize import Optimizer


class MDMin(Optimizer):
    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 dt=None):
        Optimizer.__init__(self, atoms, restart, logfile, trajectory)

        if dt is not None:
            self.dt = dt

    def initialize(self):
        self.v = None
        self.dt = 0.2

    def read(self):
        self.v, self.dt = self.load()
        
    def step(self, f):
        atoms = self.atoms

        if self.v is None:
            self.v = np.zeros((len(atoms), 3))
        else:
            self.v += 0.5 * self.dt * f
            # Correct velocities:
            vf = np.vdot(self.v, f)
            if vf < 0.0:
                self.v[:] = 0.0
            else:
                self.v[:] = f * vf / np.vdot(f, f)

        self.v += 0.5 * self.dt * f
        r = atoms.get_positions()
        atoms.set_positions(r + self.dt * self.v)
        self.dump((self.v, self.dt))
