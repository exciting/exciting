import threading

from ase.test import World
from ase.io import PickleTrajectory
from ase.neb import NEB
from ase.calculators.lj import LennardJones as Calculator
from ase.optimize import BFGS

fmax = 0.05

nimages = 3

print [a.get_potential_energy() for a in PickleTrajectory('H.traj')]
images = [PickleTrajectory('H.traj')[-1]]
for i in range(nimages):
    images.append(images[0].copy())
images[-1].positions[6, 1] = 2 - images[0].positions[6, 1]
neb = NEB(images)
neb.interpolate()

for image in images[1:]:
    image.set_calculator(Calculator())

dyn = BFGS(neb, trajectory='mep.traj')

dyn.run(fmax=fmax)

for a in neb.images:
    print a.positions[-1], a.get_potential_energy()

results = [images[2].get_potential_energy()]

def run_neb_calculation(cpu):
    images = [PickleTrajectory('H.traj')[-1]]
    for i in range(nimages):
        images.append(images[0].copy())
    images[-1].positions[6, 1] = 2 - images[0].positions[6, 1]
    neb = NEB(images, parallel=True, world=cpu)
    neb.interpolate()

    images[cpu.rank + 1].set_calculator(Calculator())

    dyn = BFGS(neb)
    dyn.run(fmax=fmax)

    if cpu.rank == 1:
        results.append(images[2].get_potential_energy())
    
w = World(nimages - 1)
ranks = [w.get_rank(r) for r in range(w.size)]
threads = [threading.Thread(target=run_neb_calculation, args=(rank,))
           for rank in ranks]
for t in threads:
    t.start()
for t in threads:
    t.join()

print results
assert abs(results[0] - results[1]) < 1e-13


