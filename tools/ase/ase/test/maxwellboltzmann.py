from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.lattice.cubic import FaceCenteredCubic

atoms = FaceCenteredCubic(size=(50,50,50), symbol="Cu", pbc=False)
print "Number of atoms:", len(atoms)
MaxwellBoltzmannDistribution(atoms, 0.1)
temp = atoms.get_kinetic_energy() / (1.5 * len(atoms))

print "Temperature", temp, " (should be 0.1)"
assert abs(temp - 0.1) < 1e-3
