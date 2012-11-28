from ase.dft import get_distribution_moment
import numpy as np

precision = 1E-8

x = np.linspace(-50., 50., 1000)
y = np.exp(-x**2 / 2.)
area, center, mom2 = get_distribution_moment(x, y, (0, 1, 2))
assert sum((abs(area - np.sqrt(2. * np.pi)), abs(center), abs(mom2 - 1.))) < precision

x = np.linspace(-1., 1., 100000)
for order in range(0, 9):
    y = x**order
    area = get_distribution_moment(x, y)
    assert abs(area - (1. - (-1.)**(order + 1)) / (order + 1.)) < precision

x = np.linspace(-50., 50., 100)
y = np.exp(-2. * (x - 7.)**2 / 10.) + np.exp(-2. * (x + 5.)**2 / 10.)
center=get_distribution_moment(x, y, 1)
assert abs(center - 1.) < precision 
