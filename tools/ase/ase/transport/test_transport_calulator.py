from ase.transport.calculators import TransportCalculator
import numpy as np

#Aux. function to write data to a text file.
def write(fname,xs,ys):
    fd = open(fname,'w')
    for x,y in zip(xs,ys):
        print >> fd, x, y
    fd.close()

H_lead = np.zeros([4,4])

# On-site energies are zero
for i in range(4):
    H_lead[i,i] = 0.0

# Nearest neighbor hopping is -1.0
for i in range(3):
    H_lead[i,i+1] = -1.0
    H_lead[i+1,i] = -1.0

# Next-nearest neighbor hopping is 0.2
for i in range(2):
    H_lead[i,i+2] = 0.2
    H_lead[i+2,i] = 0.2

H_scat = np.zeros([6,6])
# Principal layers on either side of S
H_scat[:2,:2] = H_lead[:2,:2]
H_scat[-2:,-2:] = H_lead[:2,:2]

# Scattering region
H_scat[2,2] = 0.0
H_scat[3,3] = 0.0
H_scat[2,3] = -0.8
H_scat[3,2] = -0.8

# External coupling
H_scat[1,2] = 0.2
H_scat[2,1] = 0.2
H_scat[3,4] = 0.2
H_scat[4,3] = 0.2

energies = np.arange(-3,3,0.02)
tcalc = TransportCalculator(h=H_scat,
                            h1=H_lead,
                            eta=0.02,
                            energies=energies)

T = tcalc.get_transmission()
tcalc.set(pdos=[2, 3])
pdos = tcalc.get_pdos()

tcalc.set(dos=True)
dos = tcalc.get_dos()

write('T.dat',tcalc.energies,T)
write('pdos0.dat', tcalc.energies,pdos[0])
write('pdos1.dat', tcalc.energies,pdos[1])

#subdiagonalize
h_rot, s_rot, eps, u = tcalc.subdiagonalize_bfs([2, 3], apply=True)
T_rot = tcalc.get_transmission()
dos_rot = tcalc.get_dos()
pdos_rot = tcalc.get_pdos()

write('T_rot.dat', tcalc.energies,T_rot)
write('pdos0_rot.dat', tcalc.energies, pdos_rot[0])
write('pdos1_rot.dat', tcalc.energies, pdos_rot[1])

print 'Subspace eigenvalues:', eps
assert sum(abs(eps-(-0.8, 0.8))) < 2.0e-15, 'Subdiagonalization. error'
print 'Max deviation of T after the rotation:', np.abs(T-T_rot).max()
assert max(abs(T-T_rot)) < 2.0e-15, 'Subdiagonalization. error'

#remove coupling
h_cut, s_cut = tcalc.cutcoupling_bfs([2], apply=True)
T_cut = tcalc.get_transmission()
dos_cut = tcalc.get_dos()
pdos_cut = tcalc.get_pdos()

write('T_cut.dat', tcalc.energies, T_cut)
write('pdos0_cut.dat', tcalc.energies,pdos_cut[0])
write('pdos1_cut.dat', tcalc.energies,pdos_cut[1])

