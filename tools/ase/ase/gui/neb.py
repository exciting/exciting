# -*- coding: utf-8 -*-
from math import sqrt

import numpy as np

from ase.neb import fit0

def NudgedElasticBand(images):
    N = images.repeat.prod()
    natoms = images.natoms // N

    R = images.P[:, :natoms]
    E = images.E
    F = images.F[:, :natoms]

    s, E, Sfit, Efit, lines = fit0(E, F, R)
    import pylab
    import matplotlib
    #matplotlib.use('GTK')
    pylab.ion()
    x = 2.95
    pylab.figure(figsize=(x * 2.5**0.5, x))
    pylab.plot(s, E, 'o')
    for x, y in lines:
        pylab.plot(x, y, '-g')
    pylab.plot(Sfit, Efit, 'k-')
    pylab.xlabel(u'path [Å]')
    pylab.ylabel(u'energy [eV]')
    pylab.title('Maximum: %.3f eV' % max(Efit))
    pylab.show()
