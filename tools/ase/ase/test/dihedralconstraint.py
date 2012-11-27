from ase.structure import molecule
from ase.calculators.emt import EMT
from ase.constraints import FixInternals
from ase.optimize.bfgs import BFGS

system = molecule('CH3CH2OH')
system.center(vacuum=5.0)
system.rattle(stdev=0.3)

indices = [6, 0, 1, 2]
indices2 = [6, 0, 1]
#system.set_dihedral(indices, pi/20, mask=[0,1,1,1,1,1,0,0,0])

#Angles, Bonds, Dihedrals are built up with  pairs of constraint 
#value and indices defining the constraint

angle = [system.get_angle(indices2), indices2]
dihedral = [system.get_dihedral(indices), indices] 

constraint = FixInternals(system, bonds=[], angles=[angle], dihedrals=[dihedral])

print constraint

calc = EMT()

opt = BFGS(system, trajectory='opt.traj', logfile='opt.log')

previous_angle = system.get_angle(indices2)
previous_dihedral = system.get_dihedral(indices)

print 'angle before', previous_angle
print 'dihedral before', previous_dihedral

system.set_calculator(calc)
system.set_constraint(constraint)
print '-----Optimization-----'
opt.run(fmax=0.01)

new_angle = system.get_angle(indices2)
new_dihedral = system.get_dihedral(indices)

print 'angle after', new_angle
print 'dihedral after', new_dihedral

err1 = new_angle - previous_angle
err2 = new_dihedral - previous_dihedral

print 'error in angle', repr(err1)
print 'error in dihedral', repr(err2)

assert err1 < 1e-12
assert err2 < 1e-12
