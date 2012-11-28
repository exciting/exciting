from ase import *

#from math import sqrt

#from ASE.Calculators.PairPotential import PairPotential
#from ASE.Filters.Subset import Subset
#from ASE.Filters.FixBondLength import FixBondLength
#from ASE.Dynamics.MDMin import MDMin
#from ASE.Trajectories.NetCDFTrajectory import NetCDFTrajectory
#from ASE.IO.NetCDF import ReadNetCDF

slab = read('2N.traj')
slab.set_calculator(EMT())

constraint = FixBondLength(-2, -1)
#mask=[atom.symbol != 'N' for atom in slab])

#molecule = Subset(slab, indices=[-2, -1])

# Create a trajectory for the dissociation path:
path = PickleTrajectory('association.traj', slab)
# Put the initial state in the trajectory:
path.write()

# From now on, we relax the molecule under the constraint of fixed
# bond length:
#fixed = FixBondLength(molecule, 0, 1)
#relax = MDMin(fixed, dt=0.08, fmax=0.05)

relax = QuasiNewton(slab)

d = linalg.norm(slab[-2].position - slab[-1].position)
#d = fixed.GetBondLength()
delta = 0.1
e0 = slab.get_potential_energy()
print d, 0.0
while d > 1.10:
    d -= delta
    fixed.SetBondLength(d, fix='first')
    # XXX
    raise NotImplementedError
    relax.run(fmax=.05)
    path.write()
    print d, slab.GetPotentialEnergy() - e0
