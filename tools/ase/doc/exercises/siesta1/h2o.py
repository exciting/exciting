import time
import numpy as np

from ase import Atoms, units
from ase.optimize import QuasiNewton
from ase.calculators.siesta import Siesta

# Set up a H2O molecule by specifying atomic positions
h2o = Atoms(symbols='H2O',
            positions=[( 0.776070, 0.590459, 0.00000),
                       (-0.776070, 0.590459, 0.00000),
                       (0.000000,  -0.007702,  -0.000001)],
            pbc=(1,1,1))

# Center the molecule in the cell with some vacuum around
h2o.center(vacuum=6.0)

# Select some energy-shifts for the basis orbitals
e_shifts = [0.01,0.1,0.2,0.3,0.4,0.5]

# Run the relaxation for each energy shift, and print out the
# corresponding total energy, bond length and angle


for e_s in e_shifts:
    starttime = time.time()
    calc = Siesta('h2o',meshcutoff=200.0 * units.Ry, mix=0.5, pulay=4)
    calc.set_fdf('PAO.EnergyShift', e_s * units.eV)    
    calc.set_fdf('PAO.SplitNorm', 0.15)
    calc.set_fdf('PAO.BasisSize', 'SZ')
    h2o.set_calculator(calc)
    # Make a -traj file named h2o_current_shift.traj:      
    dyn = QuasiNewton(h2o, trajectory='h2o_%s.traj' % e_s)
    dyn.run(fmax=0.02)      # Perform the relaxation      
    E = h2o.get_potential_energy()
    print                                # Make the output more readable      
    print "E_shift: %.2f" %e_s       
    print "----------------"
    print "Total Energy: %.4f" % E       # Print total energy      
    d = h2o.get_distance(0,2)
    print "Bond length: %.4f" % d        # Print bond length      
    p = h2o.positions
    d1 = p[0] - p[2]
    d2 = p[1] - p[2]
    r = np.dot(d1, d2) / (np.linalg.norm(d1) * np.linalg.norm(d2))
    angle = np.arccos(r) / pi * 180
    print "Bond angle: %.4f" % angle      # Print bond angle
    endtime = time.time()
    walltime = endtime - starttime
    print 'Wall time: %.5f' % walltime
    print                                # Make the output more readable      
