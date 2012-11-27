"""This module defines an ASE interface to Turbomole

http://www.turbomole.com/
"""
import os
import sys

import numpy as np

from ase.units import Hartree, Bohr
from ase.io.turbomole import write_turbomole
from ase.calculators.general import Calculator


class Turbomole(Calculator):
    def __init__(self, label='turbomole',
                 calculate_energy='dscf', calculate_forces='grad',
                 post_HF = False):
        self.label = label
        self.converged = False
        
        # set calculators for energy and forces
        self.calculate_energy = calculate_energy
        self.calculate_forces = calculate_forces

        # turbomole has no stress
        self.stress = np.empty((3, 3))
        
        # storage for energy and forces
        self.e_total = None
        self.forces = None
        self.updated = False
        
        # atoms must be set
        self.atoms = None
        
        # POST-HF method 
        self.post_HF  = post_HF

    def initialize(self, atoms):
        self.numbers = atoms.get_atomic_numbers().copy()
        self.species = []
        for a, Z in enumerate(self.numbers):
            self.species.append(Z)
        self.converged = False
        
    def execute(self, command):
        from subprocess import Popen, PIPE
        try:
            # the sub process gets started here
            proc = Popen([command], shell=True, stderr=PIPE)
            error = proc.communicate()[1]
            # check the error output
            if 'abnormally' in error:
                raise OSError(error)
            print 'TM command: ', command, 'successfully executed'
        except OSError, e:
            print >> sys.stderr, 'Execution failed:', e
            sys.exit(1)

    def get_potential_energy(self, atoms):
        # update atoms
        self.set_atoms(atoms)
        # if update of energy is neccessary
        if self.update_energy:
            # calculate energy
            self.execute(self.calculate_energy + ' > ASE.TM.energy.out')
            # check for convergence of dscf cycle
            if os.path.isfile('dscf_problem'):
                print 'Turbomole scf energy calculation did not converge'
                raise RuntimeError(
                    'Please run Turbomole define and come thereafter back')
            # read energy
            self.read_energy()
        else:
            print 'taking old values (E)'
        self.update_energy = False
        return self.e_total

    def get_forces(self, atoms):
        # update atoms
        self.set_atoms(atoms)
        # complete energy calculations
        if self.update_energy:
            self.get_potential_energy(atoms)
        # if update of forces is neccessary
        if self.update_forces:
            # calculate forces
            self.execute(self.calculate_forces + ' > ASE.TM.forces.out')
            # read forces
            self.read_forces()
        else:
            print 'taking old values (F)'
        self.update_forces = False
        return self.forces.copy()
    
    def get_stress(self, atoms):
        return self.stress
        
    def set_atoms(self, atoms):
        if self.atoms == atoms:
            return
        # performs an update of the atoms 
        Calculator.set_atoms(self, atoms)
        write_turbomole('coord', atoms)
        # energy and forces must be re-calculated
        self.update_energy = True
        self.update_forces = True
        
    def read_energy(self):
        """Read Energy from Turbomole energy file."""
        text = open('energy', 'r').read().lower()
        lines = iter(text.split('\n'))

        # Energy:
        for line in lines:
            if line.startswith('$end'):
                break
            elif line.startswith('$'):
                pass
            else:
                energy_tmp = float(line.split()[1])
                if self.post_HF:
                    energy_tmp += float(line.split()[4])
        # update energy units
        self.e_total = energy_tmp * Hartree

    def read_forces(self):
        """Read Forces from Turbomole gradient file."""
        file = open('gradient', 'r')
        lines = file.readlines()
        file.close()

        forces = np.array([[0, 0, 0]])
        
        nline = len(lines)
        iline = -1
        
        for i in range(nline):
            if 'cycle' in lines[i]:
                iline = i
        
        if iline < 0:
            raise RuntimeError('Please check TURBOMOLE gradients')

        # next line
        iline += len(self.atoms) + 1
        # $end line
        nline -= 1
        # read gradients
        for i in xrange(iline, nline):
            line = lines[i].replace('D', 'E')
            tmp = np.array([[float(f) for f in line.split()[0:3]]])
            forces = np.concatenate((forces, tmp))  
        # Note the '-' sign for turbomole, to get forces
        self.forces = (-np.delete(forces, np.s_[0:1], axis=0)) * Hartree / Bohr
