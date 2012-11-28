# -*- coding: utf-8 -*-
# creates:  ener.png distance.png angle.png
import os
import matplotlib
matplotlib.use('Agg')
import pylab as plt


e_s = [0.01,0.1,0.2,0.3,0.4,0.5]
E = [-463.2160, -462.9633, -462.4891, -462.0551,
     -461.5426, -461.1714]
d = [1.1131, 1.1046, 1.0960, 1.0901,
     1.0857, 1.0810]
alpha = [100.832453365, 99.568214268, 99.1486065462,
         98.873671379, 98.1726341945, 98.0535643778]

fig=plt.figure(figsize=(3, 2.5))
fig.subplots_adjust(left=.29, right=.96, top=.9, bottom=0.16)
plt.plot(e_s, E, 'o-')
plt.xlabel(u'Energy shift [eV]')
plt.ylabel(u'Energy [eV]')
plt.title('Total Energy vs Eshift')
plt.savefig('ener.png')

fig=plt.figure(figsize=(3, 2.5))
fig.subplots_adjust(left=.24, right=.96, top=.9, bottom=0.16)
plt.plot(e_s, d, 'o-')
plt.xlabel(u'Energy shift [eV]')
plt.ylabel(u'O-H distance [Å]')
limits = plt.axis('tight')
plt.title('O-H distance vs Eshift')
plt.savefig('distance.png')

fig=plt.figure(figsize=(3, 2.5))
fig.subplots_adjust(left=.26, right=.96, top=.9, bottom=0.16)
plt.plot(e_s, alpha, 'o-')
plt.xlabel(u'Energy shift [eV]')
plt.ylabel(u'H20 angle')
limits = plt.axis('tight')
plt.title('O-H distance vs Eshift')
plt.savefig('angle.png')
