"""This module defines an ASE interface to FHI-aims.

Felix Hanke hanke@liverpool.ac.uk
Jonas Bjork j.bjork@liverpool.ac.uk
"""

from general import Calculator
import os
import sys
from os.path import join, isfile, islink

import numpy as np

import ase

float_keys = [
    'charge',
    'charge_mix_param',
    'default_initial_moment',
    'hartree_convergence_parameter',
    'harmonic_length_scale',
    'ini_linear_mix_param',
    'ini_spin_mix_parma',
    'initial_moment',
    'MD_MB_init',
    'MD_time_step',
    'prec_mix_param',
    'set_vacuum_level',
    'spin_mix_param',
]

exp_keys = [
    'sc_accuracy_eev',
    'sc_accuracy_etot',
    'sc_accuracy_forces',
    'sc_accuracy_rho',
]

string_keys = [
    'communication_type',
    'density_update_method',
    'KS_method',
    'mixer',
    'output_level',
    'packed_matrix_format',
    'relax_unit_cell',
    'restart',
    'restart_read_only',
    'restart_write_only',
    'spin',
    'total_energy_method',
    'qpe_calc',
    'xc',
]

int_keys = [
    'empty_states',
    'ini_linear_mixing',
    'max_relaxation_steps',
    'multiplicity',
    'n_max_pulay',   
    'sc_iter_limit',
    'walltime',
]

bool_keys = [
    'collect_eigenvectors',
    'compute_forces',
    'compute_kinetic',
    'compute_numerical_stress',
    'distributed_spline_storage',
    'evaluate_work_function',
    'final_forces_cleaned',
    'hessian_to_restart_geometry',
    'load_balancing',
    'MD_clean_rotations',
    'MD_restart',
    'restart_relaxations',
    'squeeze_memory',
    'use_density_matrix',
    'use_dipole_correction',
    'use_local_index',
    'use_logsbt',
    'vdw_correction_hirshfeld',
]

list_keys = [
    'init_hess',
    'k_grid',
    'k_offset',
    'MD_run',
    'MD_schedule',
    'MD_segment',
    'mixer_threshold',
    'occupation_type',
    'output',
    'cube',
    'preconditioner',
    'relativistic',
    'relax_geometry',
]

input_keys = [
    'run_command',
    'run_dir',
    'species_dir',
    'cubes',
    'output_template',
    'track_output',
] 

input_parameters_default = {'run_command':None,
                            'run_dir':None,
                            'species_dir':None,
                            'cubes':None,
                            'output_template':'aims',
                            'track_output':False}

class Aims(Calculator):
    def __init__(self, **kwargs):
        self.name = 'Aims'
        self.float_params = {}
        self.exp_params = {}
        self.string_params = {}
        self.int_params = {}
        self.bool_params = {}
        self.list_params = {}
        self.input_parameters = {}
        for key in float_keys:
            self.float_params[key] = None
        for key in exp_keys:
            self.exp_params[key] = None
        for key in string_keys:
            self.string_params[key] = None
        for key in int_keys:
            self.int_params[key] = None
        for key in bool_keys:
            self.bool_params[key] = None
        for key in list_keys:
            self.list_params[key] = None
        for key in input_keys:
            self.input_parameters[key] = input_parameters_default[key]
        if os.environ.has_key('AIMS_SPECIES_DIR'):
            self.input_parameters['species_dir'] = os.environ['AIMS_SPECIES_DIR']
        if os.environ.has_key('AIMS_COMMAND'):
            self.input_parameters['run_command'] = os.environ['AIMS_COMMAND']

        self.positions = None
        self.atoms = None
        self.run_counts = 0
        self.set(**kwargs)

    def set(self, **kwargs):
        if 'control' in kwargs:
            fname = kwargs['control']
            from ase.io.aims import read_aims_calculator
            file = open(fname, 'r')
            calc_temp = None
            while True:
                line = file.readline()
                if "List of parameters used to initialize the calculator:" in line:
                   file.readline()
                   calc_temp = read_aims_calculator(file)
                   break
            if calc_temp is not None:
                self.float_params      = calc_temp.float_params
                self.exp_params        = calc_temp.exp_params
                self.string_params     = calc_temp.string_params
                self.int_params        = calc_temp.int_params
                self.bool_params       = calc_temp.bool_params
                self.list_params       = calc_temp.list_params
                self.input_parameters  = calc_temp.input_parameters
            else:
                raise TypeError("Control file to be imported can not be read by ASE = " + fname)
        for key in kwargs:
            if self.float_params.has_key(key):
                self.float_params[key] = kwargs[key]
            elif self.exp_params.has_key(key):
                self.exp_params[key] = kwargs[key]
            elif self.string_params.has_key(key):
                self.string_params[key] = kwargs[key]
            elif self.int_params.has_key(key):
                self.int_params[key] = kwargs[key]
            elif self.bool_params.has_key(key):
                self.bool_params[key] = kwargs[key]
            elif self.list_params.has_key(key):
                self.list_params[key] = kwargs[key]
            elif self.input_parameters.has_key(key):
                self.input_parameters[key] = kwargs[key]
            elif key is not 'control':
                raise TypeError('Parameter not defined: ' + key)

    def update(self, atoms):
        if self.calculation_required(atoms,[]):
            self.calculate(atoms)

    def calculation_required(self, atoms,quantities):
        if (self.positions is None or
            (self.atoms != atoms) or
            (self.atoms != self.old_atoms) or 
            (self.float_params != self.old_float_params) or
            (self.exp_params != self.old_exp_params) or
            (self.string_params != self.old_string_params) or
            (self.int_params != self.old_int_params) or
            (self.bool_params != self.old_bool_params) or
            (self.list_params != self.old_list_params) or
            (self.input_parameters != self.old_input_parameters)):
            return True
        else:
            return False

    def calculate(self, atoms):
        """Generate necessary files in the working directory.
        
        If the directory does not exist it will be created.

        """
        positions = atoms.get_positions()
        have_lattice_vectors = atoms.get_pbc().any()        
        have_k_grid = self.list_params['k_grid']
        if have_lattice_vectors and not have_k_grid:
            raise RuntimeError("Found lattice vectors but no k-grid!")
        if not have_lattice_vectors and have_k_grid:
            raise RuntimeError("Found k-grid but no lattice vectors!")
        from ase.io.aims import write_aims
        write_aims('geometry.in', atoms) 
        self.write_control()
        self.write_species()
        self.run()
        self.converged = self.read_convergence()
        if not self.converged:
            os.system("tail -20 "+self.out)
            raise RuntimeError("FHI-aims did not converge!\n"+
                               "The last lines of output are printed above "+
                               "and should give an indication why.")

        self.set_results(atoms)

    def set_results(self,atoms):
        self.read(atoms)
        self.old_float_params = self.float_params.copy()
        self.old_exp_params = self.exp_params.copy()
        self.old_string_params = self.string_params.copy()
        self.old_int_params = self.int_params.copy()
        self.old_input_parameters = self.input_parameters.copy()
        self.old_bool_params = self.bool_params.copy()
        self.old_list_params = self.list_params.copy()
        self.old_atoms = self.atoms.copy()

    def run(self):
        if self.input_parameters['track_output']:
            self.out = self.input_parameters['output_template']+str(self.run_counts)+'.out'
            self.run_counts += 1
        else:
            self.out = self.input_parameters['output_template']+'.out'            
        if self.input_parameters['run_command']:
            aims_command = self.input_parameters['run_command'] 
        elif os.environ.has_key('AIMS_COMMAND'):
            aims_command = os.environ['AIMS_COMMAND']
        else:
            raise RuntimeError("No specification for running FHI-aims. Aborting!")
        aims_command = aims_command + ' >> ' 
        if self.input_parameters['run_dir']:
            aims_command = aims_command + self.input_parameters['run_dir'] + '/'
        aims_command = aims_command + self.out
        self.write_parameters('#',self.out)
        exitcode = os.system(aims_command)
        if exitcode != 0:
            raise RuntimeError('FHI-aims exited with exit code: %d.  ' % exitcode)
        if self.input_parameters['cubes'] and self.input_parameters['track_output']:
            self.input_parameters['cubes'].move_to_base_name(self.input_parameters['output_template']+str(self.run_counts-1))

    def write_parameters(self,prefix,filename):
        output = open(filename,'w')
        output.write(prefix+'=======================================================\n')
        output.write(prefix+'FHI-aims file: '+filename+'\n')
        output.write(prefix+'Created using the Atomic Simulation Environment (ASE)\n'+prefix+'\n')
        output.write(prefix+'List of parameters used to initialize the calculator:\n')
        output.write(prefix+'=======================================================\n')
        for key, val in self.float_params.items():
            if val is not None:
                output.write('%-35s%5.6f\n' % (key, val))        
        for key, val in self.exp_params.items():
            if val is not None:
                output.write('%-35s%5.2e\n' % (key, val))
        for key, val in self.string_params.items():
            if val is not None:
                output.write('%-35s%s\n' % (key, val))
        for key, val in self.int_params.items():
            if val is not None:
                output.write('%-35s%d\n' % (key, val))
        for key, val in self.bool_params.items():
            if val is not None:
                if key == 'vdw_correction_hirshfeld' and val:
                    output.write('%-35s\n' % (key))
                elif val:
                    output.write('%-35s.true.\n' % (key))
                elif key != 'vdw_correction_hirshfeld':
                    output.write('%-35s.false.\n' % (key))
        for key, val in self.list_params.items():
            if val is not None:
                if key == 'output':
                    if not isinstance(val,(list,tuple)): 
                        val = [val]
                    for output_type in val:
                        output.write('%-35s%s\n' % (key,str(output_type)))
                else:
                    output.write('%-35s' % (key))
                    if isinstance(val,str): 
                        output.write(val)
                    else:
                        for sub_value in val:
                            output.write(str(sub_value)+' ')
                    output.write('\n')
        for key, val in self.input_parameters.items():
            if key is  'cubes':
                if val:
                    val.write(output)
            elif val and val != input_parameters_default[key]:
                output.write(prefix+'%-34s%s\n' % (key,val))
        output.write(prefix+'=======================================================\n\n')
        output.close()

    def __repr__(self):
        items =  self.float_params.items()+self.exp_params.items() \
                +self.string_params.items()+self.int_params.items()
        rep = 'Aims('
        for key, val in items:
            if val is not None:
                rep += key+' '+str(val)+', '
        for key, val in self.bool_params.items():
            if val is not None:
                if key == 'vdw_correction_hirshfeld' and val:
                    rep += key + ', '
                elif val:
                    rep += key+' .true., '
                elif key != 'vdw_correction_hirshfeld':
                    rep += key+' .false., '
        for key, val in self.list_params.items():
            if val is not None:
                if key == 'output':
                    if not isinstance(val,(list,tuple)): 
                        val = [val]
                    for output_type in val:
                        rep += key+' '+str(output_type)+', '
                else:
                    rep += key
                    if isinstance(val,str): 
                        rep += ' '+val
                    else:
                        for sub_value in val:
                            rep += ' '+str(sub_value)
                    rep += ', '
        for key, val in self.input_parameters.items():
            if val and val != input_parameters_default[key]:
                rep += key+' '+val+', '
        return rep[:-2]+')'
        
    def write_control(self, file = 'control.in'):
        """Writes the control.in file."""
        self.write_parameters('#',file)

    def write_species(self, file = 'control.in'):
        from ase.data import atomic_numbers

        if not self.input_parameters['species_dir']:
            raise RuntimeError('Missing species directory, THIS MUST BE SPECIFIED!')

        control = open(file, 'a')
        species_path = self.input_parameters['species_dir']
        symbols = self.atoms.get_chemical_symbols()
        symbols2 = []
        for n, symbol in enumerate(symbols):
            if symbol not in symbols2:
                symbols2.append(symbol)
        for symbol in symbols2:
            fd = join(species_path, '%02i_%s_default' % (atomic_numbers[symbol], symbol))
            for line in open(fd, 'r'):
                control.write(line)
        control.close()

    def get_dipole_moment(self, atoms):
        if self.list_params['output'] is None or 'dipole' not in self.list_params['output']:
            raise RuntimeError('output=[\'dipole\'] has to be set.')
        elif atoms.get_pbc().any():
            raise RuntimeError('FHI-aims does not allow this for systems with periodic boundary conditions.')
        self.update(atoms)
        return self.dipole

    def read_dipole(self):
        """Method that reads the electric dipole moment from the output file."""

        dipolemoment=np.zeros([1,3])
        for line in open(self.out, 'r'):
            if line.rfind('Total dipole moment [eAng]') > -1:
                dipolemoment=np.array([float(f) for f in line.split()[6:10]])
        return dipolemoment

    def read_energy(self, all=None):
        for line in open(self.out, 'r'):
            if line.rfind('Total energy corrected') > -1:
                E0 = float(line.split()[5])
            elif line.rfind('Total energy uncorrected') > -1:
                F = float(line.split()[5])
        energy_free, energy_zero = F, E0
        return [energy_free, energy_zero]

    def read_forces(self, atoms, all=False):
        """Method that reads forces from the output file.

        If 'all' is switched on, the forces for all ionic steps
        in the output file will be returned, in other case only the
        forces for the last ionic configuration are returned."""
        lines = open(self.out, 'r').readlines()
        forces = np.zeros([len(atoms), 3])
        for n, line in enumerate(lines):
            if line.rfind('Total atomic forces') > -1:
                for iatom in range(len(atoms)):
                    data = lines[n+iatom+1].split()
                    for iforce in range(3):
                        forces[iatom, iforce] = float(data[2+iforce])
        return forces

    def read_stress(self):
        lines = open(self.out, 'r').readlines()
        stress = None
        for n, line in enumerate(lines):
            if line.rfind('Calculation of numerical stress completed') > -1:
                stress = []
                for i in [n+8,n+9,n+10]:
                    data = lines[i].split()
                    stress += [float(data[2]),float(data[3]),float(data[4])]
        # rearrange in 6-component form and return
        if stress is not None:
            return np.array([stress[0], stress[4], stress[8], stress[5], stress[2], stress[1]])
        else:
            return

    def get_stress(self, atoms):
        self.update(atoms)
        return self.stress

# methods that should be quickly implemented some time, haven't had time yet:
    def read_fermi(self):
        """Method that reads Fermi energy from output file"""
        return

    def read_magnetic_moment(self):
        return

    def read_convergence(self):
        converged = False
        lines = open(self.out, 'r').readlines()
        for n, line in enumerate(lines):
            if line.rfind('Have a nice day') > -1:
                converged = True
        return converged

    def read_eigenvalues(self, kpt=0, spin=0):
        return 

class AimsCube:
    """ object to ensure the output of cube files, can be attached to Aims object"""
    def __init__(self,origin=(0,0,0),
                 edges=[(0.1,0.0,0.0),(0.0,0.1,0.0),(0.0,0.0,0.1)],
                 points=(50,50,50),plots=None):
        """ parameters: 
        origin, edges, points = same as in the FHI-aims output
        plots: what to print, same names as in FHI-aims """

        self.name   = 'AimsCube'
        self.origin = origin
        self.edges  = edges
        self.points = points
        self.plots  = plots
         
    def ncubes(self):
        """returns the number of cube files to output """
        if self.plots:
            number = len(self.plots)
        else:
            number = 0
        return number

    def set(self,**kwargs):
        """ set any of the parameters ... """
        # NOT IMPLEMENTED AT THE MOMENT!

    def move_to_base_name(self,basename):
        """ when output tracking is on or the base namem is not standard,
        this routine will rename add the base to the cube file output for 
        easier tracking """
        for plot in self.plots:
            found = False
            cube = plot.split()
            if cube[0] == 'total_density' or cube[0] == 'spin_density' or cube[0] == 'delta_density':
                found = True
                old_name = cube[0]+'.cube'
                new_name = basename+'.'+old_name
            if cube[0] == 'eigenstate' or cube[0] == 'eigenstate_density':
                found = True
                state = int(cube[1])
                s_state = cube[1]
                for i in [10,100,1000,10000]:
                    if state < i:
                        s_state = '0'+s_state
                old_name = cube[0]+'_'+s_state+'_spin_1.cube'
                new_name = basename+'.'+old_name
            if found:
                os.system("mv "+old_name+" "+new_name)

    def add_plot(self,name):
        """ in case you forgot one ... """
        plots += [name]

    def write(self,file):
        """ write the necessary output to the already opened control.in """
        file.write('output cube '+self.plots[0]+'\n')
        file.write('   cube origin ')
        for ival in self.origin:
            file.write(str(ival)+' ')
        file.write('\n')
        for i in range(3):
            file.write('   cube edge '+str(self.points[i])+' ')
            for ival in self.edges[i]:
                file.write(str(ival)+' ')
            file.write('\n')
        if self.ncubes() > 1:
            for i in range(self.ncubes()-1):
                file.write('output cube '+self.plots[i+1]+'\n')

                    
                
