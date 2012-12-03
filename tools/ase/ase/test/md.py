from ase import Atoms
from ase.calculators.emt import EMT
from ase.md import VelocityVerlet
from ase.io import PickleTrajectory

a = 3.6
b = a / 2
fcc = Atoms('Cu', positions=[(0, 0, 0)],
            cell=[(0, b, b), (b, 0, b), (b, b, 0)],
            pbc=1)
fcc *= (2, 1, 1)
fcc.set_calculator(EMT())
fcc.set_momenta([(0.9, 0.0, 0.0), (-0.9, 0, 0)])
md = VelocityVerlet(fcc, dt=0.1)
def f():
    print fcc.get_potential_energy(), fcc.get_total_energy()
md.attach(f)
md.attach(PickleTrajectory('Cu2.traj', 'w', fcc).write, interval=3)
md.run(steps=20)
fcc2 = PickleTrajectory('Cu2.traj', 'r')[-1]

