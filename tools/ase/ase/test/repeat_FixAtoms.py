from ase.structure import molecule
from ase.constraints import FixAtoms

N = 2

atoms = molecule('CO2')
atoms.set_cell((15.,15.,15.))

print('indices method')
atomsi = atoms.copy()
atomsi.set_constraint(FixAtoms(indices=[0,]))
atomsi = atomsi.repeat((N,1,1))

atomsiref = atoms.copy().repeat((N,1,1))
atomsiref.set_constraint(FixAtoms(indices=[0, N + 1]))

lcatomsi = list(atomsi.constraints[0].index)
lcatomsiref = list(atomsiref.constraints[0].index)

assert lcatomsi == lcatomsiref

print('mask method')
atomsm = atoms.copy()
atomsm.set_constraint(FixAtoms(mask=[True, False, False]))
atomsm = atomsm.repeat((N,1,1))

atomsmref = atoms.copy().repeat((N,1,1))
atomsmref.set_constraint(FixAtoms(mask=[True, False, False] * N))

lcatomsm = list(atomsm.constraints[0].index)
lcatomsmref = list(atomsmref.constraints[0].index)

assert lcatomsm == lcatomsmref

# http://stackoverflow.com/questions/3873361/finding-multiple-occurrences-of-a-string-within-a-string-in-python

lcatomsm2i = [n for (n, e) in enumerate(lcatomsm) if e == True]

assert lcatomsm2i == lcatomsi
