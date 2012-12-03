# -*- coding: utf-8 -*-
from math import sqrt

import numpy as np

from ase.units import kJ
from ase.utils.eos import EquationOfState


def BulkModulus(images):
    v = np.array([abs(np.linalg.det(A)) for A in images.A])
    #import matplotlib.pyplot as plt
    import pylab as plt
    plt.ion()
    EquationOfState(v, images.E).plot()

"""
    fit = np.poly1d(np.polyfit(v**-(1.0 / 3), images.E, 3))
    fit1 = np.polyder(fit, 1)
    fit2 = np.polyder(fit1, 1)
    for t in np.roots(fit1):
        if t > 0 and fit2(t) > 0:
            break
    v0 = t**-3
    e0 = fit(t)
    B = t**5 * fit2(t) / 9 / kJ * 1.0e24  # Gpa

    import pylab
    import matplotlib
    #matplotlib.use('GTK')

    pylab.ion()
    x = 3.95
    pylab.figure(figsize=(x * 2.5**0.5, x))

    pylab.plot(v, images.E, 'o')

    x = np.linspace(min(v), max(v), 100)
    pylab.plot(x, fit(x**-(1.0 / 3)), '-r')
    pylab.xlabel(u'volume [Å^3]')
    pylab.ylabel(u'energy [eV]')
    pylab.title(u'E: %.3f eV, V: %.3f Å^3, B: %.3f GPa' % (e0, v0, B))
    pylab.show()
"""
