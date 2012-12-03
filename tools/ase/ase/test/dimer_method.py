from ase.lattice.surface import fcc100, add_adsorbate
from ase.constraints import FixAtoms
from ase.calculators.emt import EMT
from ase.dimer import DimerControl, MinModeAtoms, MinModeTranslate

# Set up a small "slab" with an adatoms
atoms = fcc100('Pt', size = (2, 2, 1), vacuum = 10.0)
add_adsorbate(atoms, 'Pt', 1.611, 'hollow')

# Freeze the "slab"
mask = [atom.tag > 0 for atom in atoms]
atoms.set_constraint(FixAtoms(mask = mask))

# Calculate using EMT
atoms.set_calculator(EMT())
relaxed_energy = atoms.get_potential_energy()

# Set up the dimer
d_control = DimerControl(initial_eigenmode_method = 'displacement', \
                         displacement_method = 'vector', logfile = None, \
                         mask = [0, 0, 0, 0, 1])
d_atoms = MinModeAtoms(atoms, d_control)

# Dispalce the atoms
displacement_vector = [[0.0]*3]*5
displacement_vector[-1][1] = -0.1
d_atoms.displace(displacement_vector = displacement_vector)

# Converge to a saddle point
dim_rlx = MinModeTranslate(d_atoms, trajectory = 'dimer_method.traj', \
                           logfile = None)
dim_rlx.run(fmax = 0.001)

# Test the results
tolerance = 1e-3
assert(d_atoms.get_barrier_energy() - 1.03733136918 < tolerance)
assert(abs(d_atoms.get_curvature() + 0.900467048707) < tolerance)
assert(d_atoms.get_eigenmode()[-1][1] < -0.99)
assert(abs(d_atoms.get_positions()[-1][1]) < tolerance)
