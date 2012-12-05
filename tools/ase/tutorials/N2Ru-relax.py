from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton
from ase.constraints import FixAtoms

a = 2.70
c = 1.59 * a
h = 1.85
d = 1.10

slab = Atoms('2Cu', [(0., 0., 0.), (1/3., 1/3., -0.5*c)], 
             tags=(0, 1),
             pbc=(1, 1, 0))
slab.set_cell([(a, 0, 0),
               (a / 2, 3**0.5 * a / 2, 0),
               (0, 0, 1)])
slab = slab.repeat((4, 4, 1))
slab.set_calculator(EMT())
mask = [a.tag == 1 for a in slab]
slab.set_constraint(FixAtoms(mask=mask))
dyn = QuasiNewton(slab)
dyn.run(fmax=0.05)

e_slab = slab.get_potential_energy()
x = slab.positions[0, 2] / (c / 2) * 100
print 'Relaxation of the top layer: %f %%' % x

molecule = Atoms('2N', positions=[(0., 0., h),
                                  (0., 0., h + d)])
molecule.set_calculator(EMT())
e_N2 = molecule.get_potential_energy()
slab.extend(molecule)

dyn = QuasiNewton(slab)
dyn.run(fmax=0.05)

print 'Adsorption energy:', e_slab + e_N2 - slab.get_potential_energy()
