from ase.lattice.surface import fcc100, add_adsorbate
from ase.constraints import FixAtoms, FixedPlane
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton

# 2x2-Al(001) surface with 3 layers and an
# Au atom adsorbed in a hollow site:
slab = fcc100('Al', size=(2, 2, 3))
add_adsorbate(slab, 'Au', 1.7, 'hollow')
slab.center(axis=2, vacuum=4.0)

# Make sure the structure is correct:
#view(slab)

# Fix second and third layers:
mask = [atom.tag > 1 for atom in slab]
#print mask
fixlayers = FixAtoms(mask=mask)

# Constrain the last atom (Au atom) to move only in the yz-plane:
plane = FixedPlane(-1, (1, 0, 0))

slab.set_constraint([fixlayers, plane])

# Use EMT potential:
slab.set_calculator(EMT())

# Initial state:
for i in range(5):
    qn = QuasiNewton(slab, trajectory='mep%d.traj' % i)
    qn.run(fmax=0.05)
    slab[-1].x += slab.get_cell()[0, 0] / 8
