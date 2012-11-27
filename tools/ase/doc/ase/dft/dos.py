# creates: dos.png

import numpy as np

from ase.dft import DOS

class MyCalc:
    def get_eigenvalues(self, kpt=0, spin=0):
        return np.random.uniform(-5.0, 2.0, 90)
    def get_k_point_weights(self):
        return [1.0]
    def get_number_of_spins(self):
        return 1
    def get_fermi_level(self):
        return 0.0

calc = MyCalc()
dos = DOS(calc, width=0.2)
d = dos.get_dos()
e = dos.get_energies()

import matplotlib
matplotlib.use('Agg')
import pylab as plt
plt.figure(figsize=(5, 4))
plt.plot(e, d)
plt.xlabel('energy [eV]')
plt.ylabel('DOS')
plt.savefig('dos.png')


