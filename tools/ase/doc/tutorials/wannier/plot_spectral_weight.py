from ase.dft import Wannier
from gpaw import restart

atoms, calc = restart('benzene.gpw', txt=None)
wan = Wannier(nwannier=18, calc=calc, fixedstates=15, file='wan18.pickle')

import pylab as pl
weight_n = pl.sum(abs(wan.V_knw[0])**2, 1)
N = len(weight_n)
F = wan.fixedstates_k[0]
pl.figure(1, figsize=(12, 4))
pl.bar(range(1, N+1), weight_n, width=0.65, bottom=0,
        color='k', edgecolor='k', linewidth=None,
       align='center', orientation='vertical')
pl.plot([F+.5, F+.5], [0, 1], 'k--')
pl.axis(xmin=.32, xmax=N+1.33, ymin=0, ymax=1)
pl.xlabel('Eigenstate')
pl.ylabel('Projection of wannier functions')
pl.savefig('spectral_weight.png')
pl.show()
