from ase import *

# XXX This cannot be converted to ase 3, since we lack ANEB
raise NotImplementedError

#from ASE.Calculators.PairPotential import PairPotential
#from ASE.Filters.Subset import Subset
#from ASE.Dynamics.AdaptiveNudgedElasticBand import AdaptiveNudgedElasticBand
#from ASE.IO.NetCDF import ReadNetCDF
import os

# check that N2.nc and 2N.nc exists otherwise run
# N2Ru-Dissociation1.py
if not (os.path.exists('N2.nc') and os.path.exists('N2.nc')):
    os.system('python N2Ru-Dissociation1.py')

initial = read('N2.traj')
final = read('2N.traj')

configs = [initial.copy() for i in range(4)] + [final]

constraint = FixAtoms(mask=[atom.symbol != 'N' for atom in initial])
for config in configs:
    config.set_calculator(EMT())
    config.set_constraint(constraint)
    
band = NEB(configs)
band.interpolate()

# Create a quickmin object:
relax = QuasiNewton(band)



#configs0 = [initial]
#for i in range(3):
#    configs0.append(initial.Copy())
#configs0.append(final)
#
#configs = []
#for config in configs0:
#    config.SetCalculator(PairPotential())
#    configs.append(Subset(config, indices=[-2, -1]))


# setup the Adaptive Nudged Elastic Band dynamics
aneb = AdaptiveNudgedElasticBand(configs,prefix='aneb',linear_interpolation=True,maxlevels=1)

aneb.Converge()

# test the anebfit tool
os.system('anebfit aneb')
os.system('xmgrace aneb_-1.agr&')

# clean up
# os.system('rm aneb*conf*')

