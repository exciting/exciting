#!/usr/bin/python
"""This simple demo calculates the total energy of CO molecules
using once LDA and once PBE as xc-functional. Obviously
some parts in this scripts are longer than necessary, but are shown
to demonstrate some more features."""

import ase
import ase.calculators.castep, ase.io.castep

calc = ase.calculators.castep.Castep()
directory = 'CASTEP_ASE_DEMO'

# include interface settings in .param file
calc._export_settings = True

# reuse the same directory
calc._directory = directory
calc._rename_existing_dir = False
calc._label = 'CO_LDA'

# necessary for tasks with changing positions
# such as GeometryOptimization or MolecularDynamics
calc._set_atoms = True

# Param settings
calc.param.xc_functional = 'LDA'
calc.param.cut_off_energy = 400
# Prevent CASTEP from writing *wvfn* files
calc.param.num_dump_cycles = 0

# Cell settings
calc.cell.kpoint_mp_grid = '1 1 1'
calc.cell.fix_com = False
calc.cell.fix_all_cell = True

# Set and clear and reset settings (just for shows)
calc.param.task = 'SinglePoint'
# Reset to CASTEP default
calc.param.task.clear()

# all of the following are identical
calc.param.task = 'GeometryOptimization'
calc.task = 'GeometryOptimization'
calc.TASK = 'GeometryOptimization'
calc.Task = 'GeometryOptimization'


# Prepare atoms
mol = ase.atoms.Atoms('CO', [[0, 0, 0], [0, 0, 1.2]], cell=[10, 10, 10])
mol.set_calculator(calc)

# Check for correct input
if calc.dryrun_ok():
    print('%s : %s ' % (mol.calc._label, mol.get_potential_energy()))
else:
    print("Found error in input")
    print(calc._error)


# Read all settings from previous calculation
mol = ase.io.castep.read_seed('%s/CO_LDA' % directory)

# Use the OTF pseudo-potential we have just generated
mol.calc.set_pspot('OTF')

# Change some settings
mol.calc.param.xc_functional = 'PBE'
# don't forget to set an appropriate label
mol.calc._label = 'CO_PBE'
# Recalculate the potential energy
print('%s : %s ' % (mol.calc._label, mol.get_potential_energy()))
