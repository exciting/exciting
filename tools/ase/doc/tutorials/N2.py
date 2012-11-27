from ase import Atoms
from ase.calculators.emt import EMT

atom = Atoms('N', calculator=EMT())
e_atom = atom.get_potential_energy()

d = 1.1
molecule = Atoms('2N', [(0., 0., 0.), (0., 0., d)])
molecule.set_calculator(EMT())
e_molecule = molecule.get_potential_energy()

e_atomization = e_molecule - 2 * e_atom

print 'Nitrogen atom energy: %5.2f eV' %e_atom
print 'Nitrogen molecule energy: %5.2f eV' %e_molecule
print 'Atomization energy: %5.2f eV' %-e_atomization
