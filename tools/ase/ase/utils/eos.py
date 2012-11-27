# -*- coding: utf-8 -*-
from math import sqrt

import numpy as np

from ase.units import kJ

class EquationOfState:
    """Fit equation of state for bulk systems.

    The following equation is used::

                          2      3        -1/3
      E(V) = c + c t + c t  + c t ,  t = V
              0   1     2      3

    Use::

       eos = EquationOfState(volumes, energies)
       v0, e0, B = eos.fit()
       eos.plot()

    """
    def __init__(self, volumes, energies):
        self.v = np.array(volumes)
        self.e = np.array(energies)
        self.v0 = None

    def fit(self):
        """Calculate volume, energy, and bulk modulus.

        Returns the optimal volume, the minumum energy, and the bulk
        modulus.  Notice that the ASE units for the bulk modulus is
        eV/Angstrom^3 - to get the value in GPa, do this::

          v0, e0, B = eos.fit()
          print B / kJ * 1.0e24, 'GPa'
          
        """
        
        fit0 = np.poly1d(np.polyfit(self.v**-(1.0 / 3), self.e, 3))
        fit1 = np.polyder(fit0, 1)
        fit2 = np.polyder(fit1, 1)

        self.v0 = None
        for t in np.roots(fit1):
            if t > 0 and fit2(t) > 0:
                self.v0 = t**-3
                break

        if self.v0 is None:
            raise ValueError('No minimum!')
        
        self.e0 = fit0(t)
        self.B = t**5 * fit2(t) / 9
        self.fit0 = fit0
        
        return self.v0, self.e0, self.B

    def plot(self, filename=None, show=None):
        """Plot fitted energy curve.

        Uses Matplotlib to plot the energy curve.  Use *show=True* to
        show the figure and *filename='abc.png'* or
        *filename='abc.eps'* to save the figure to a file."""
        
        #import matplotlib.pyplot as plt
        import pylab as plt

        if self.v0 is None:
            self.fit()
            
        if filename is None and show is None:
            show = True

        x = 3.95
        f = plt.figure(figsize=(x * 2.5**0.5, x))
        f.subplots_adjust(left=0.12, right=0.9, top=0.9, bottom=0.15)
        plt.plot(self.v, self.e, 'o')
        x = np.linspace(min(self.v), max(self.v), 100)
        plt.plot(x, self.fit0(x**-(1.0 / 3)), '-r')
        plt.xlabel(u'volume [Å^3]')
        plt.ylabel(u'energy [eV]')
        plt.title(u'E: %.3f eV, V: %.3f Å^3, B: %.3f GPa' %
                  (self.e0, self.v0, self.B / kJ * 1.0e24))

        if show:
            plt.show()
        if filename is not None:
            f.savefig(filename)

        return f
