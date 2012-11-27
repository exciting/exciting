from ase.data.molecules import molecule
from gpaw import GPAW

atoms = molecule('C6H6')
atoms.center(vacuum=3.5)

calc = GPAW(h=.21, xc='PBE', txt='benzene.txt', nbands=18)
atoms.set_calculator(calc)
atoms.get_potential_energy()

calc.set(fixdensity=True, txt='benzene-harris.txt',
         nbands=40, eigensolver='cg', convergence={'bands': 35})
atoms.get_potential_energy()

calc.write('benzene.gpw', mode='all')
