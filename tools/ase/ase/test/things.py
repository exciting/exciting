from ase.dft import monkhorst_pack

assert [0, 0, 0] in  monkhorst_pack((1, 3, 5)).tolist()
assert [0, 0, 0] not in  monkhorst_pack((1, 3, 6)).tolist()
assert len(monkhorst_pack((3, 4, 6))) == 3 * 4 * 6

from ase.units import Hartree, Bohr, kJ, mol, kcal, kB, fs
print Hartree, Bohr, kJ/mol, kcal/mol, kB*300, fs, 1/fs

from ase.lattice import bulk
hcp = bulk('X', 'hcp', a=1) * (2, 2, 1)
assert abs(hcp.get_distance(0, 3, mic=True) - 1) < 1e-12
assert abs(hcp.get_distance(0, 4, mic=True) - 1) < 1e-12
assert abs(hcp.get_distance(2, 5, mic=True) - 1) < 1e-12
