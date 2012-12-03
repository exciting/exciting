from math import pi, sqrt

import numpy as np


class DOS:
    def __init__(self, calc, width=0.1, window=None, npts=201):
        """Electronic Density Of States object.

        calc: calculator object
            Any ASE compliant calculator object.
        width: float
            Width of guassian smearing.
        window: tuple of two float
            Use ``window=(emin, emax)``.  If not specified, a window
            big enough to hold all the eigenvalues will be used.
        npts: int
            Number of points.

        """
        
        self.npts = npts
        self.width = width
        self.w_k = calc.get_k_point_weights()
        self.nspins = calc.get_number_of_spins()
        self.e_skn = np.array([[calc.get_eigenvalues(kpt=k, spin=s)
                                for k in range(len(self.w_k))]
                               for s in range(self.nspins)])
        self.e_skn -= calc.get_fermi_level()

        if window is None:
            emin = self.e_skn.min() - 5 * self.width
            emax = self.e_skn.max() + 5 * self.width
        else:
            emin, emax = window

        self.energies = np.linspace(emin, emax, npts)

    def get_energies(self):
        """Return the array of energies used to sample the DOS. The energies are reported relative to the Fermi level."""
        return self.energies

    def delta(self, energy):
        """Return a delta-function centered at 'energy'."""
        x = -((self.energies - energy) / self.width)**2
        return np.exp(x) / (sqrt(pi) * self.width)

    def get_dos(self, spin=None):
        """Get array of DOS values.

        The *spin* argument can be 0 or 1 (spin up or down) - if not
        specified, the total DOS is returned.
        """
        
        if spin is None:
            if self.nspins == 2:
                # Spin-polarized calculation, but no spin specified -
                # return the total DOS:
                return self.get_dos(spin=0) + self.get_dos(spin=1)
            else:
                spin = 0
        
        dos = np.zeros(self.npts)
        for w, e_n in zip(self.w_k, self.e_skn[spin]):
            for e in e_n:
                dos += w * self.delta(e)
        return dos
