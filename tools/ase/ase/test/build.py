import numpy as np
from ase import Atoms, Atom

a = Atoms([Atom('Cu')])
a.positions[:] += 1.0
print a.get_positions(), a.positions
a=a+a
a+=a
a.append(Atom('C'))
a += Atoms([])
a += Atom('H', magmom=1)
print a.get_initial_magnetic_moments()
print a[0].number
print a[[0,1]].get_atomic_numbers()
print a[np.array([1,1,0,0,1], bool)].get_atomic_numbers()
print a[::2].get_atomic_numbers()
print a.get_chemical_symbols()
del a[2]
print a.get_chemical_symbols()
del a[-2:]
print a.get_chemical_symbols()
