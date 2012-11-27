"""This module defines an ASE interface to NWchem

http://www.nwchem-sw.org/
"""
import os
import sys

import numpy as np

from ase.units import Hartree, Bohr
from ase.io.nwchem import write_nwchem
from ase.calculators.general import Calculator

class KPoint:
    def __init__(self, s):
        self.s = s
        self.eps_n = []
        self.f_n = []

class NWchem(Calculator):
    def __init__(self,
                 label='nwchem',
                 task='energy',
                 # Warning: nwchem centers atoms by default
                 # see ase-developers/2012-March/001356.html
                 geometry='nocenter',
                 xc='LDA',
                 convergence = {'energy'  : None,
                                'density' : None,
                                'gradient': None,
                                'lshift': None, # set to 0.0 for nolevelshifting
                                },
                 smear=0.001*Hartree, # smear must be specified to appear in out!
                 grid=None,
                 tolerances=None,
                 cgmin=False,
                 maxiter = 120,
                 basis='3-21G',
                 basispar=None,
                 ecp=None,
                 so=None,
                 charge=None,
                 multiplicity=None,
                 spinorbit=False,
                 kpts=None,
                 dftcontrol='', # additional dft control string
                 control='', # additional outside of dft block control string
                 ):
        """Construct NWchem-calculator object.

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.nw, label.out, ...).
            Default is 'nwchem'.
        xc: str
            Exchange-correlation functional. LDA and PBE are predefined,
            use nchem names instead.
        basis: str
            Basis set
        maxiter: int
            Maximal number of iteratations in self-consistent field convergence.
            """

        self.label = label
        self.task = task
        self.geometry = geometry
        self.xc = xc
        self.convergence = convergence
        self.smear = round(smear/Hartree, 4)
        self.grid = grid
        self.tolerances = tolerances
        self.cgmin = cgmin
        self.maxiter = maxiter
        self.basis = basis
        if basispar is not None:
            self.basispar = 'basis ' + basispar
        else:
            self.basispar = 'basis'
        self.ecp = ecp
        self.so = so
        self.charge = charge
        self.multiplicity = multiplicity
        self.spinorbit = spinorbit
        self.kpts = kpts
        self.dftcontrol = dftcontrol
        self.control = control

        # does nwchem have stress ???
        self.stress = np.zeros((3, 3))

        # atoms must be set
        self.atoms = None

        self.converged = False

    def execute(self, command):
        from subprocess import Popen, PIPE
        try:
            # the sub process gets started here
            proc = Popen([command], shell=True, stderr=PIPE)
            error = proc.communicate()[1]
            if error:
                raise OSError(error + '\ncheck ' + self.output)
        except OSError, e:
            print >> sys.stderr, 'Execution failed:', e
            sys.exit(1)

    def run(self):
        """Method which explicitely runs Nwchem."""

        command = os.environ.get('NWCHEM_COMMAND', 'nwchem')
        # The label may contain escapable characters!
        self.execute(command + ' "' + \
                     self.label + '.nw" > "' + self.output + '"')

    def get_forces(self, atoms):
        self.get_potential_energy(atoms)
        return self.forces

    def get_ibz_k_points(self):
        return np.array([0., 0., 0.])

    def get_electronic_temperature(self):
        return self.electronic_temperature

    def read_smear(self):
        smear = None
        for line in open(self.label + '.out'): # find last one
            if line.find('Smearing applied:') != -1:
                smear = float(line.split(':')[1].split()[0].strip().lower().replace('d', 'e'))
        return smear

    def get_number_of_bands(self):
        return self.nvector

    def read_number_of_bands(self):
        nvector = 0
        for line in open(self.label + '.out'):
            if line.find('Vector ') != -1: # count all printed vectors
                nvector += 1
        if not nvector:
            nvector = None
        return nvector

    def get_number_of_electrons(self):
        return self.nelect

    def read_number_of_electrons(self):
        nelect = None
        for line in open(self.label + '.out'): # find last one
            if line.find('of electrons') != -1:
                nelect = float(line.split(':')[1].strip())
        return nelect

    def get_number_of_iterations(self):
        return self.niter

    def read_number_of_iterations(self):
        niter = 0
        for line in open(self.label + '.out'):
            if line.find('d= ') != -1: # count all iterations
                niter += 1
        if not niter:
            niter = None
        return niter

    def get_magnetic_moment(self, atoms):
        return self.magnetic_moment

    def read_magnetic_moment(self):
        magmom = None
        for line in open(self.label + '.out'):
            if line.find('Spin multiplicity') != -1: # last one
                magmom = float(line.split(':')[-1].strip()) - 1
        return magmom

    def get_magnetic_moments(self, atoms):
        # local magnetic moments are not available in nwchem
        # so set the total magnetic moment on the atom no. 0 and fill with 0.0
        magmoms = [0.0 for a in range(len(atoms))]
        magmoms[0] = self.get_magnetic_moment(atoms)
        return np.array(magmoms)

    def get_dipole_moment(self, atoms=None):
        return self.dipole

    def read_dipole_moment(self):
        dipolemoment=[]
        for line in open(self.label + '.out'):
            for component in [
                '1   1 0 0',
                '1   0 1 0',
                '1   0 0 1'
                ]:
                if line.find(component) != -1:
                    value = float(line.split(component)[1].split()[0])  # total dipole component
                    value = value * Bohr
                    dipolemoment.append(value)
        if len(dipolemoment) == 0:
            dipolemoment = None
        return dipolemoment

    def get_potential_energy(self, atoms):
        # update atoms
        self.set_atoms(atoms)
        # if update of energy is neccessary
        if self.energy is None or self.forces is None:
            # write input file
            f = open(self.label + '.nw', 'w')
            if self.charge is not None:
                f.write('charge ' + str(self.charge) + '\n')
            write_nwchem(f, atoms, self.geometry)

            def format_basis_set(string, tag=self.basispar):
                formatted = tag + '\n'
                lines = string.split('\n')
                if len(lines) > 1:
                    formatted += string
                else:
                    formatted += '  * library '  + string + '\n'
                return formatted + 'end\n'
            basis = format_basis_set(self.basis)
            if self.ecp is not None:
                basis += format_basis_set(self.ecp, 'ecp')
            if self.so is not None:
                basis += format_basis_set(self.so, 'so')
            f.write(basis)

            if self.xc == 'RHF':
                task = 'scf'
            else:
                if self.spinorbit:
                    task = 'sodft'
                else:
                    task = 'dft'
                nwchem_xc_map = {
                    'LDA' : 'slater pw91lda',
                    'PBE' : 'xpbe96 cpbe96',
                    }
                if self.xc in nwchem_xc_map:
                    xc = nwchem_xc_map[self.xc]
                else:
                    xc = self.xc
                f.write('\n' + task + '\n')
                f.write('  mult ' + str(self.multiplicity) + '\n')
                f.write('  xc ' + xc + '\n')
                f.write('  iterations ' + str(self.maxiter) + '\n')
                for key in self.convergence:
                    if key == 'lshift':
                        if self.convergence[key] is not None:
                            if not (self.convergence[key] > 0.0):
                                f.write('  convergence nolevelshifting\n')
                            else:
                                f.write('  convergence ' + key + ' ' +
                                        str(self.convergence[key]/Hartree) + '\n')
                    else:
                        if self.convergence[key] is not None:
                            f.write('  convergence ' + key + ' ' +
                                    str(self.convergence[key]) + '\n')
                if self.smear is not None:
                    f.write('  smear ' + str(self.smear) + '\n')
                if self.grid is not None:
                    f.write('  grid ' + str(self.grid) + '\n')
                if self.tolerances is not None:
                    f.write('  tolerances ' + str(self.tolerances) + '\n')
                if self.cgmin:
                    f.write('  cgmin\n')
                if self.dftcontrol:
                    f.write(self.dftcontrol + '\n')
                f.write('end\n')

            if self.control:
                f.write(self.control + '\n')

#            f.write('\ntask ' + task + ' gradient\n')
            f.write('\ntask ' + task + ' ' + self.task + '\n')
            f.close()

            # calculate energy
            self.output = self.label + '.out'
            self.run()
            # read output
            self.read_energy()
            if self.task.find('gradient') > -1:
                self.read_forces()
            self.niter = self.read_number_of_iterations()
            self.nelect = self.read_number_of_electrons()
            self.nvector = self.read_number_of_bands()
            self.magnetic_moment = self.read_magnetic_moment()
            self.electronic_temperature = self.read_smear() * Hartree
            self.dipole = self.read_dipole_moment()
        else:
            print 'taking old values (E)'

        return self.energy * Hartree

    def read_energy(self):
        """Read Energy from nwchem output file."""
        text = open(self.output, 'r').read()
        lines = iter(text.split('\n'))

        # Energy:
        for line in lines:
            estring = 'Total '
            if self.xc == 'RHF':
                estring += 'SCF'
            else:
                estring += 'DFT'
            estring += ' energy'
            if line.find(estring) >=0:
                energy = float(line.split()[4])
                break
        self.energy = energy

        # Eigenstates
        spin = -1
        kpts = []
        for line in lines:
            if line.find('Molecular Orbital Analysis') >= 0:
                spin += 1
                kpts.append(KPoint(spin))
            if spin >= 0:
                if line.find('Vector') >= 0:
                    line = line.lower().replace('d', 'e')
                    line = line.replace('=', ' ')
                    word = line.split()
                    kpts[spin].f_n.append(float(word[3]))
                    kpts[spin].eps_n.append(float(word[5]))
        self.kpts = kpts

    def read_forces(self):
        """Read Forces from nwchem output file."""
        file = open(self.output, 'r')
        lines = file.readlines()
        file.close()

        for i, line in enumerate(lines):
            if line.find('ENERGY GRADIENTS') >=0:
                gradients = []
                for j in range(i + 4, i + 4 + len(self.atoms)):
                    word = lines[j].split()
                    gradients.append([float(word[k]) for k in range(5,8)])
        self.forces =  - np.array(gradients) * Hartree / Bohr

    def get_eigenvalues(self, kpt=0, spin=0):
        """Return eigenvalue array."""
        return np.array(self.kpts[spin].eps_n) * Hartree

    def get_occupation_numbers(self, kpt=0, spin=0):
        """Return occupation number array."""
        return self.kpts[spin].f_n

    def get_number_of_spins(self):
        """Return the number of spins in the calculation.

        Spin-paired calculations: 1, spin-polarized calculation: 2."""
        return len(self.kpts)

    def get_spin_polarized(self):
        """Is it a spin-polarized calculation?"""
        return len(self.kpts) == 2

    def set_atoms(self, atoms):
        if self.atoms == atoms:
            return

        self.atoms = atoms.copy()
        self.energy = None
        self.forces = None

        if self.multiplicity is None:
            # obtain multiplicity from magnetic momenta
            multiplicity = 1 + atoms.get_initial_magnetic_moments().sum()
            self.multiplicity = int(multiplicity)
            if self.multiplicity != multiplicity:
                raise RuntimeError('Noninteger multiplicity not possible.\n' +
                                   'Check initial magnetic moments.')
        else:
            if self.multiplicity != int(self.multiplicity):
                raise RuntimeError('Noninteger multiplicity not possible.')

    def update(self, atoms):
        self.set_atoms(atoms)
