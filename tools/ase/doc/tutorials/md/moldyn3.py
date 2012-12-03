"Demonstrates molecular dynamics with constant energy."

from ase.calculators.emt import EMT
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.langevin import Langevin
from ase.io.trajectory import PickleTrajectory
from ase import units

from asap3 import EMT   # Way too slow with ase.EMT !
size = 10

T = 1500 # Kelvin

# Set up a crystal
atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]], symbol="Cu",
                          size=(size,size,size), pbc=False)

# Describe the interatomic interactions with the Effective Medium Theory
atoms.set_calculator(EMT())

# We want to run MD with constant energy using the Langevin algorithm
# with a time step of 5 fs, the temperature T and the friction
# coefficient to 0.02 atomic units.
dyn = Langevin(atoms, 5*units.fs, T*units.kB, 0.002)

#Function to print the potential, kinetic and total energy.
def printenergy(a=atoms):    #store a reference to atoms in the definition.
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print ("Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  Etot = %.3feV" %
           (epot, ekin, ekin/(1.5*units.kB), epot+ekin))
dyn.attach(printenergy, interval=50)

#We also want to save the positions of all atoms after every 100th time step.
traj = PickleTrajectory("moldyn3.traj", 'w', atoms)
dyn.attach(traj.write, interval=50)

# Now run the dynamics
printenergy()
dyn.run(5000)

