"Demonstrates molecular dynamics for isolated particles."

from ase.calculators.emt import EMT
from ase.cluster.cubic import FaceCenteredCubic
from ase.optimize import QuasiNewton
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, \
                                        Stationary, ZeroRotation
from ase.md.verlet import VelocityVerlet
from ase import units

# Use Asap for a huge performance increase if it is installed
useAsap = True

if useAsap:
    from asap3 import EMT
    size = 4
else:
    size = 2

# Set up a nanoparticle
atoms = FaceCenteredCubic('Cu', surfaces=[[1,0,0],[1,1,0],[1,1,1]],
                          layers=(size,size,size), vacuum=4)

# Describe the interatomic interactions with the Effective Medium Theory
atoms.set_calculator(EMT())

# Do a quick relaxation of the cluster
qn = QuasiNewton(atoms)
qn.run(0.001, 10)

# Set the momenta corresponding to T=1200K
MaxwellBoltzmannDistribution(atoms, 1200*units.kB)
Stationary(atoms) # zero linear momentum
ZeroRotation(atoms) # zero angular momentum

# We want to run MD using the VelocityVerlet algorithm.
dyn = VelocityVerlet(atoms, 5*units.fs, trajectory='moldyn4.traj') # save trajectory.

#Function to print the potential, kinetic and total energy.
def printenergy(a=atoms):    #store a reference to atoms in the definition.
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print ("Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  Etot = %.3feV" %
           (epot, ekin, ekin/(1.5*units.kB), epot+ekin))
dyn.attach(printenergy, interval=10)

# Now run the dynamics
printenergy()
dyn.run(2000)
