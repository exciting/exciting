import numpy as np

class DOS:

    def __init__(self, energy: np.ndarray, dos: np.ndarray):
        self.energy = energy
        self.dos = dos
