# creates: Al_phonon.png Al_mode.gif Al_mode.pdf

from ase.lattice import bulk
from ase.calculators.emt import EMT
from ase.dft.kpoints import ibz_points, get_bandpath
from ase.phonons import Phonons

# Setup crystal and EMT calculator
atoms = bulk('Al', a=4.05)
calc = EMT()

# Phonon calculator
N = 6
ph = Phonons(atoms, calc, supercell=(N, N, N))
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True)

# High-symmetry points in the Brillouin zone
points = ibz_points['fcc']
G = points['Gamma']
X = points['X']
W = points['W']
K = points['K']
L = points['L']
U = points['U']

point_names = ['$\Gamma$', 'X', 'U', 'L', '$\Gamma$', 'K']
path = [G, X, U, L, G, K]
path_kc, q, Q = get_bandpath(path, atoms.cell, 100)
omega_kn = 1000 * ph.band_structure(path_kc)

# DOS
omega_e, dos_e = ph.dos(kpts=(50, 50, 50), npts=5000, delta=1e-4)
omega_e *= 1000

# Plot phonon dispersion
import matplotlib
matplotlib.use('Agg')
import pylab as plt

plt.figure(1, (8, 6))
plt.axes([.1, .07, .67, .85])
for n in range(len(omega_kn[0])):
    omega_n = omega_kn[:, n]
    plt.plot(q, omega_n, 'k-', lw=2)

plt.xticks(Q, point_names, fontsize=18)
plt.yticks(fontsize=18)
plt.xlim(q[0], q[-1])
plt.ylim(0, 35)
plt.ylabel("Frequency ($\mathrm{meV}$)", fontsize=22)
plt.grid('on')

plt.axes([.8, .07, .17, .85])
plt.fill_between(dos_e, omega_e, y2=0, color='lightgrey', edgecolor='k', lw=2)
plt.ylim(0, 35)
plt.xticks([], [])
plt.yticks([], [])
plt.xlabel("DOS", fontsize=18)
plt.savefig('Al_phonon.png')

# Write modes for specific q-vector to trajectory files
ph.write_modes([l/2 for l in L], branches=[2], repeat=(8, 8, 8), kT=3e-4,
               center=True)

# Generate png animation
from subprocess import call
from ase.io import PickleTrajectory, write

trajfile = 'phonon.mode.2.traj'
trajectory = PickleTrajectory(trajfile, 'r')

for i, atoms in enumerate(trajectory):
    write('picture%02i.png' %i, atoms, show_unit_cell=2,
          rotation='-36x,26.5y,-25z')
    # Flatten images for better quality
    call(['convert', '-flatten', 'picture%02i.png' %i, 'picture%02i.png' %i])

# Make static pdf image for pdflatex
call(['convert', 'picture00.png', 'Al_mode.pdf'])

# Concatenate to gif animation
call(['convert', '-delay', '5', '-loop', '0', '-dispose', 'Previous', 'picture*.png',
      'Al_mode.gif'])



