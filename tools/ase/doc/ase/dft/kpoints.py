# creates: cc.png
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pylab as plt
from ase.dft.kpoints import cc162_1x1
B = [(1, 0, 0), (-0.5, 3**0.5 / 2, 0), (0, 0, 1)]
k = np.dot(cc162_1x1, B)
plt.figure(figsize=(5, 4))
plt.plot(k[:, 0], k[:, 1], 'o')
plt.savefig('cc.png')


