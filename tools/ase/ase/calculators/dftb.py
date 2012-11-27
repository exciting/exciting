"""This module defines an ASE interface to DftbPlus

http://http://www.dftb-plus.info//
http://www.dftb.org/

-Markus Kaukonen markus.kaukonen@iki.fi
"""


import numpy as np
import os

#from ase.data import chemical_symbols
from ase.units import Hartree, Bohr


class Dftb:
    """Class for doing DFTB+ calculations.
    """
    def __init__(self, label='dftb', write_dftb=False,
                 charge=0.0, include_dispersion=False,
                 do_spin_polarized=False, 
                 unpaired_electrons=0.0,
                 fermi_temperature=0.0, scc=False,
                 kpt_n11=1, kpt_n12=0, kpt_n13=0,
                 kpt_n21=0, kpt_n22=1, kpt_n23=0,
                 kpt_n31=0, kpt_n32=0, kpt_n33=1,
                 kpt_s1=0.0, kpt_s2=0.0, kpt_s3=0.0,
                 md_dftb=False, dftb_nsteps=0, md_restart_frequency=100,
                 md_dftb_write_vel=True, md_dftb_init_t=300.0,
                 md_dftb_berendsen_nvt=False,md_dftb_berendsen_t=300,
                 md_dftb_keep_stationary=True,
                 md_dftb_dt=0.1,md_dftb_berendsen_strength=1.0e-4):
        """Construct DFTB-calculator object.


        For example:
        calc = Dftb(label='dftb',write_dftb=True,include_dispersion=True )

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.txt, ...).
            Default is 'dftb'.
        write_dftb: boolean
            True: a minimal input file (name of which is always 'dftb_in.hsd')
            is written based on values given here.
            False: input file for dftb+ is not written. User must have
            generated file 'dftb_in.hsd' in the working directory.
            Use write_dftb=False to use your own 'dftb_in.hsd'-file.
        charge: float
            Total charge of the system.
        include_dispersion: boolean
            True: Default dispersion parameters are written in the 
            file 'dftb_in.hsd' (requires that also write_dftb_input_file==True)
            False: dispersion parameters are not written here.
        do_spin_polarized: boolean
            True: Spin polarized calculation
            (requires that also write_dftb_input_file==True)
            False: Spin unpolarized calculation
        unpaired_electrons: float
            Number of spin unpaired electrons in the system.
            Relevant only if do_spin_polarized==True
        fermi_temperature: float
            Fermi temperature for electrons.
        scc: boolean
            True: Do charge self consistent dftb+
            False: No SCC, charges on atoms are not iterated
        kpt_nxx: int(s)
            The 9 integers for K-point generation scheme SupercellFolding
        kpt_s1, kpt_s2, kpt_s3: float(s)
            Shifts for the K-point generation scheme SupercellFolding
        md_dftb: boolean
            True: MD is run by DFTB
            False: MD is run by ase (or A CG run is done by ase or DFTB)
            (requires that also write_dftb_input_file==True)
        dftb_nsteps: int
            Number of simulation steps (either in CG or MD)
            run by DFTB+ program.
            If dftb_nsteps==0 it is assumed that minimisation/dynamics is run
            by ase (after DFTB+ positions and velocities are read).
            If dftb_nsteps>0 minimisation/dynamics is run by DFTB+.
            (requires that also write_dftb_input_file==True)
        md_restart_frequency: int
            How often coordinates are written.
            (requires that md_dftb==True, ie MD is run by DFTB+ program)
            (requires that also write_dftb_input_file==True)
        md_dftb_write_vel: bool
            True:
                velocities are written in dftb_in.hsd file to be read by
                dftb (in this case md_dftb_init_t is ignored)
            False:
                velocities are not written in dftb_in.hsd
                (requires that md_dftb==True, ie MD is run by DFTB+ program)
                (requires that also write_dftb_input_file==True)    
        md_dftb_init_t: float
            Initial temperature is set randomly to md_dftb_init_t
            This is dummy if md_dftb_write_vel==True
            (requires that md_dftb==True, ie MD is run by DFTB+ program)
            (requires that also write_dftb_input_file==True)            
        md_dftb_berendsen_nvt: boolean
            True: Berendsen NVT dynamics is run by DFTB+ (for dftb_nsteps)
            False: no Berendsen NVT MD
            (requires that also write_dftb_input_file==True)
        md_dftb_berendsen_t: float
            Temperature for DFTB+ Berendsen thermostate
            (requires that also write_dftb_input_file==True)
        md_dftb_keep_stationary: boolean
            True: In DFTB+ MD the translational movement
            of the system is removed
            False: Center of mass movement is not removed in DFTB+ MD
            (requires that also write_dftb_input_file==True)
        md_dftb_dt: float
            Timestep in [fs] for DFTB+ Berendsen thermostate.
            (requires that also write_dftb_input_file==True)
        md_dftb_berendsen_strength:float
            Coupling strength in DFTB+ Berendsen thermostate
            (requires that also write_dftb_input_file==True)
            
        Input file for DFTB+ file is 'dftb_in.hsd'. Keywords in it are
        written here or read from an existing file. The atom positions
        in file 'dftb_in.hsd' are updated during ASE geometry
        optimization.
        """


        if not(write_dftb):
            if os.path.isfile('dftb_in.hsd'):
                myf = open('dftb_in.hsd')
            else:
                print 'Input file for DFTB+ dftb_in.hsd is missing'
                raise RuntimeError, \
                    'Provide it or set write_dftb=True '
            #lines = f.readlines()
            myf.close()

        self.label = label
        self.write_dftb = write_dftb
        self.charge = charge
        self.include_dispersion = include_dispersion
        self.do_spin_polarized = do_spin_polarized
        self.unpaired_electrons = unpaired_electrons
        #if (do_spin_polarized):
        #    print 'Sorry, generation of file "dftb_in.hsd" with spin '
        #    print 'polarization is not inplemented for DFTB+'
        #    print 'Set write_dftb=False and'
        #    raise RuntimeError, \
        #        'Generate file "dftb_in.hsd by hand"'
        
        self.etotal = 0.0
        self.cell = None
        self.fermi_temperature = fermi_temperature
        self.scc = scc
        self.kpt_n11 = kpt_n11
        self.kpt_n12 = kpt_n12
        self.kpt_n13 = kpt_n13
        self.kpt_n21 = kpt_n21
        self.kpt_n22 = kpt_n22
        self.kpt_n23 = kpt_n23
        self.kpt_n31 = kpt_n31
        self.kpt_n32 = kpt_n32
        self.kpt_n33 = kpt_n33
        self.kpt_s1 = kpt_s1
        self.kpt_s2 = kpt_s2
        self.kpt_s3 = kpt_s3
        self.md_dftb = md_dftb
        self.dftb_nsteps = dftb_nsteps
        self.md_restart_frequency = md_restart_frequency
        self.md_dftb_write_vel = md_dftb_write_vel
        self.md_dftb_init_t = md_dftb_init_t
        self.md_dftb_berendsen_nvt = md_dftb_berendsen_nvt        
        self.md_dftb_berendsen_t = md_dftb_berendsen_t
        self.md_dftb_keep_stationary = md_dftb_keep_stationary
        self.md_dftb_dt = md_dftb_dt
        self.md_dftb_berendsen_strength = md_dftb_berendsen_strength
        self.converged = False

        #dftb has no stress
        self.stress = np.empty((3, 3))
        

    def update(self, atoms):
        """Energy and forces are calculated when atoms have moved
        by calling self.calculate
        """
        if (not self.converged or
            len(self.typenumber) != len(atoms)):
            self.initialize(atoms)
            self.calculate(atoms)
        elif ((self.positions != atoms.get_positions()).any() or
              (self.pbc != atoms.get_pbc()).any() or
              (self.cell != atoms.get_cell()).any()):
            self.calculate(atoms)

    def initialize(self, atoms):
        """ Set arrays of species, number of atoms in species and 
            max angular momentum for each species.
            Also write dftb_in.hsd file if desired
        """
        #from ase.io.dftb import read_dftb

        atomtypes = atoms.get_chemical_symbols()
        self.allspecies = []
        self.typenumber = []
        self.max_angular_momentum = []
        for species in (atomtypes):
            if species not in self.allspecies:
                self.allspecies.append(species)
        for species in (atomtypes):
            myindex = 1 + self.allspecies.index(species)
            self.typenumber.append(myindex)
        for i in self.allspecies:
            if i == 'H':
                self.max_angular_momentum.append('s')
            elif i in ['C','N','O']:
                self.max_angular_momentum.append('p')
            elif i in ['Si','S','Fe','Ni']:
                self.max_angular_momentum.append('d')
            else:
                print 'anglular momentum is not imlemented in ASE-DFTB+'
                print 'for species '+i
                raise RuntimeError('Use option write_dftb=False')
        self.converged = False
        
        #write DFTB input file if desired 
        self.positions = atoms.get_positions().copy()
        if self.write_dftb:
            self.write_dftb_input_file(atoms)


    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.etotal

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()
    
    def get_stress(self, atoms):
        self.update(atoms)
        return self.stress.copy()

    def calculate(self, atoms):
        """Total DFTB energy is calculated (to file 'energy'
        also forces are calculated (to file 'gradient')
        In case of md_dftb==True (md is run by DFTB+)
        the positions and velocities are updated after the DFTB+ md run.
        """
        
        self.positions = atoms.get_positions().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()

        if self.write_dftb:
            self.write_dftb_input_file(atoms)
        else:
            #write current coordinates to file 'dftb_in.hsd' for DFTB+
            self.change_atom_positions_dftb(atoms)

        #print "OK1"
        #DFTB energy&forces calculation (or possibly full MD if md_dftb==True)
        if (os.environ.has_key('DFTB_COMMAND') and
            (os.environ.has_key('DFTB_PREFIX'))):
            dftb = os.environ['DFTB_COMMAND']
            exitcode = os.system(dftb+ '> dftb.output' )
        elif not(os.environ.has_key('DFTB_COMMAND')):
            raise RuntimeError('Please set DFTB_COMMAND environment variable')
        elif not(os.environ.has_key('DFTB_PREFIX')):
            print 'Path for DFTB+ slater koster files is missing'    
            raise RuntimeError('Please set DFTB_PREFIX environment variable')
        if exitcode != 0:
            raise RuntimeError('Dftb exited with exit code: %d.  ' % exitcode)

        #read positions and velocities in case MD/CG was run by DFTB
        if (self.md_dftb):
            #print "READINGG POSITIONS"
            self.read_positions(atoms)
            self.read_velocities(atoms)
            #print "setting vel after DFTB+",atoms.get_velocities().copy()
        else:
            # in case of dynamics/minimisation is done by ASE
            # read energies&forces
            self.read_energy()
            #DFTB atomic forces calculation, to be read in file detailed.out
            #os.system(self.dftb_program_forces +'> output.forces.dummy')
            self.read_forces(atoms)
            #in case of cg minimisation with more than 0 steps 
            #we read the DFTB+ minimised geometry
            if (self.dftb_nsteps > 0):
                self.read_positions(atoms)
            
        self.converged = True

        
    def read_energy(self):
        """Read Energy from DFTB energy file."""
        text = open('dftb.output', 'r').read().lower()
        lines = iter(text.split('\n'))

        # Energy:
        for line in lines:
            if 'total energy' in line:
                energy_tmp = float(line.split()[2])
                #print 'energy_tmp', energy_tmp
        self.etotal = energy_tmp * Hartree


    def read_forces(self, atoms):
        """Read Forces from DFTB+ detailed.out output file"""

        myf = open('detailed.out','r')
        line = myf.readline()
        line = myf.readline()
        tmpforces = np.array([[0, 0, 0]])
        while line:
            if 'Total Forces' in line:
                for dummy in atoms:
                    line = myf.readline()
                    line2 = str.replace(line, 'D', 'E')
                    tmp = np.array\
                        ([[float(f) for f in line2.split()[0:3]]])
                    tmpforces = np.concatenate((tmpforces, tmp))  
            line = myf.readline()
            

        self.forces = (np.delete(tmpforces, np.s_[0:1], axis=0))*Hartree/Bohr

    def read_positions(self, atoms):
        """Read positions from LABEL.xyz output file"""

        myf = open(self.label+'.xyz','r')

        # take last coordinates of the LABEL.xyz file
        lines = myf.readlines()
        linesreversed = []
        for linereverse in reversed(lines):
            if ('MD iter' in linereverse):
                break
            else:
                linesreversed.append(linereverse)
        tmpcoord = np.array([[0, 0, 0]])
        for lineok in reversed(linesreversed):
            tmp = np.array\
                  ([[float(f) for f in lineok.split()[1:4]]])
            tmpcoord = np.concatenate((tmpcoord, tmp))  
        atoms.set_positions(np.delete(tmpcoord, np.s_[0:1], axis=0))

    def read_velocities(self, atoms):
        """Read velocities from LABEL.xyz output file"""

        myf = open(self.label+'.xyz','r')

        # take last coordinates of the LABEL.xyz file
        lines = myf.readlines()
        linesreversed = []
        for linereverse in reversed(lines):
            if ('MD iter' in linereverse):
                break
            else:
                linesreversed.append(linereverse)
        tmpvel = np.array([[0, 0, 0]])
        for lineok in reversed(linesreversed):
            tmp = np.array\
                  ([[float(f) for f in lineok.split()[4:7]]])
            tmpvel = np.concatenate((tmpvel, tmp))  
        tmpvel = tmpvel / 98.22693531550318
        atoms.set_velocities(np.delete(tmpvel, np.s_[0:1], axis=0))
        print "velocities after dftb", atoms.get_velocities().copy()


    def read(self):
        """Dummy stress for dftb"""
        self.stress = np.empty((3, 3))

    def write_dftb_input_file(self, atoms):
        """Write input parameters to DFTB+ input file 'dftb_in.hsd'."""
        import math

        myf = open('dftb_in.hsd', 'w')
        # geometry
        myf.write('Geometry = {\n')
        myf.write('TypeNames = {')
        for i in self.allspecies:
            myf.write(' "'+i+'"')
        myf.write(' }\n')
        myf.write('TypesAndCoordinates [Angstrom] = {\n')
        self.positions = atoms.get_positions().copy()
        for i, pos in zip(self.typenumber, self.positions):
            myf.write('%6d ' % (i))
            myf.write('%20.14f %20.14f %20.14f' %  tuple(pos))
            myf.write('\n')
        myf.write(' }\n')

        #is it periodic and when is write also lattice vectors
        periodic = atoms.get_pbc().any()
        if periodic:
            myf.write('Periodic = Yes\n')
        else:
            myf.write('Periodic = No\n')
        if periodic:
            cell = atoms.get_cell().copy()
            myf.write('LatticeVectors [Angstrom] = {\n')
            for latvec in cell:
                myf.write('%20.14f %20.14f %20.14f \n' %  tuple(latvec))
            myf.write('  }\n')

        #end of geometry session
        myf.write('}\n')

        #print self.dftb_nsteps
        if self.md_dftb:
            #Md run by DFTB
            myf.write('\n')
            myf.write('Driver = VelocityVerlet{\n')
            myf.write('  Steps = '+str(self.dftb_nsteps)+' \n')
            if self.md_dftb_keep_stationary:
                myf.write('  KeepStationary = Yes \n')
            else:
                myf.write('  KeepStationary = No \n')
            myf.write('  TimeStep [fs] = '+str(self.md_dftb_dt)+' \n')

            if self.md_dftb_berendsen_nvt:
                myf.write('  Thermostat = Berendsen{ \n')
                myf.write('  Temperature [K] = '+
                    str(self.md_dftb_berendsen_t)+ '\n')
                myf.write('    CouplingStrength = '+
                    str(self.md_dftb_berendsen_strength)+ ' \n')
                myf.write('  }\n')
            else:
                myf.write('  Thermostat = None{ \n')
                if not self.md_dftb_write_vel:
                    myf.write('    InitialTemperature [K] = '+
                             str(self.md_dftb_init_t)+ '\n')
                myf.write('  }\n')

                
            myf.write('  OutputPrefix = '+self.label+ '\n')
            myf.write('  MDRestartFrequency = '+
                     str(self.md_restart_frequency)+ '\n')
            if self.md_dftb_write_vel:
                myf.write('  Velocities [AA/ps]={ \n')
                tmpvelocities = atoms.get_velocities().copy()
                #from ase units to dftb+ units Ang/ps
                tmpvelocities = tmpvelocities*98.22693531550318
                for vel in tmpvelocities:
                    if math.isnan(vel[0]) or \
                            math.isnan(vel[1]) or \
                            math.isnan(vel[2]):
                        myf.write('  %20.14f %20.14f %20.14f' %  (0, 0, 0))
                    else:
                        myf.write('  %20.14f %20.14f %20.14f' %  tuple(vel))
                    myf.write('\n')
                myf.write('  }\n')
            myf.write('}\n')
        else:
            #zero step CG relaxation to get forces and energies
            myf.write('\n') 
            myf.write('Driver = ConjugateGradient {\n')
            myf.write('MovedAtoms = Range { 1 -1 }\n')
            myf.write('  MaxForceComponent = 1.0e-4\n')
            myf.write('  MaxSteps = '+str(self.dftb_nsteps)+' \n')
            myf.write('  OutputPrefix = '+self.label+ '\n')
            myf.write('}\n')

        #Hamiltonian
        myf.write('\n') 
        myf.write('Hamiltonian = DFTB { # DFTB Hamiltonian\n')
        if (self.scc):
            myf.write('  SCC = Yes')
            myf.write(' # Use self consistent charges\n')               
            myf.write('  SCCTolerance = 1.0e-5')
            myf.write(' # Tolerance for charge consistence\n')            
            myf.write('  MaxSCCIterations = 50')
            myf.write(' # Nr. of maximal SCC iterations\n')          
            myf.write('  Mixer = Broyden {') 
            myf.write(' # Broyden mixer for charge mixing\n')          
            myf.write('    MixingParameter = 0.2')  
            myf.write(' # Mixing parameter\n')
            myf.write('  }\n')
        else:
            myf.write('  SCC = No # NO self consistent charges\n')
        myf.write('  SlaterKosterFiles = Type2FileNames {')
        myf.write(' # File names from two atom type names\n')
        sk_prefix = os.environ['DFTB_PREFIX']
        myf.write('    Prefix = "'+sk_prefix+'"')
        myf.write(' # Path as prefix\n')
        myf.write('    Separator = "-"')
        myf.write(' # Dash between type names\n')
        myf.write('    Suffix = ".skf"')
        myf.write(' # Suffix after second type name\n')
        myf.write('  }\n')
        myf.write('  MaxAngularMomentum = {')
        myf.write(' # Maximal l-value of the various species\n')
        for i, j in zip(self.allspecies, self.max_angular_momentum):
            myf.write('   '+i+' = "'+j+'"\n')
        myf.write('  }\n')
        myf.write('  Charge = ')
        myf.write('%10.6f' % (self.charge))
        myf.write(' # System neutral\n')
        if self.do_spin_polarized:
            myf.write('  SpinPolarisation = Colinear {\n') 
            myf.write('  UnpairedElectrons = ' + \
                          str(self.unpaired_electrons)+'\n')
            myf.write('  } \n')
            myf.write('  SpinConstants = {\n') 
            for i in self.allspecies:
                if i == 'H':
                    myf.write('   H={\n') 
                    myf.write('    # Wss\n') 
                    myf.write('    -0.072\n') 
                    myf.write('    }\n')
                elif i == 'C': 
                    myf.write('   C={\n') 
                    myf.write('    # Wss Wsp Wps Wpp\n') 
                    myf.write('    -0.031 -0.025 -0.025 -0.023\n') 
                    myf.write('    }\n') 
                elif i == 'N': 
                    myf.write('   N={\n') 
                    myf.write('    # Wss Wsp Wps Wpp\n')
                    myf.write('    -0.033 -0.027 -0.027 -0.026\n') 
                    myf.write('     }\n') 
                elif i == 'O':
                    myf.write('   O={\n') 
                    myf.write('    # Wss Wsp Wps Wpp\n') 
                    myf.write('    -0.035 -0.030 -0.030 -0.028\n') 
                    myf.write('     }\n') 
                elif (i == 'Si' or i == 'SI'):
                    myf.write('   Si={\n') 
                    myf.write('    # Wss Wsp Wsd Wps Wpp Wpd Wds Wdp Wdd\n')
                    myf.write('    -0.020 -0.015 0.000 -0.015 \
-0.014 0.000 0.002 0.002 -0.032\n')
                    myf.write('    }\n')
                elif (i == 'S'):
                    myf.write('   S={\n') 
                    myf.write('    # Wss Wsp Wsd Wps Wpp Wpd Wds Wdp Wdd\n') 
                    myf.write('    -0.021 -0.017 0.000 -0.017 \
-0.016 0.000 0.000 0.000 -0.080\n')
                    myf.write('    }\n')
                elif (i == 'Fe' or i == 'FE'):
                    myf.write('   Fe={\n') 
                    myf.write('    # Wss Wsp Wsd Wps Wpp Wpd Wds Wdp Wdd\n') 
                    myf.write('    -0.016 -0.012 -0.003 -0.012 \
-0.029 -0.001 -0.003 -0.001 -0.015\n')
                    myf.write('    }\n')
                elif (i == 'Ni' or i == 'NI'):
                    myf.write('   Ni={\n') 
                    myf.write('    # Wss Wsp Wsd Wps Wpp Wpd Wds Wdp Wdd\n') 
                    myf.write('    -0.016 -0.012 -0.003 -0.012 \
-0.022 -0.001 -0.003 -0.001 -0.018\n')
                    myf.write('    }\n')
                else:
                    print 'Missing spin polarisation parameters for species'+i
                    raise RuntimeError, \
                        'Run spin unpolarised calculation'
            myf.write('   }\n') 
        else:
            #myf.write('  SpinPolarisation = {}')
            myf.write(' # No spin polarisation\n')
        myf.write('  Filling = Fermi {\n')
        myf.write('    Temperature [Kelvin] = ')
        myf.write('%10.6f\n' % (self.fermi_temperature))
        myf.write('  }\n')
        if periodic:
            myf.write('  KPointsAndWeights = SupercellFolding { \n')
            myf.write('%6d %6d %6d \n' 
                     % (self.kpt_n11, self.kpt_n12, self.kpt_n13))
            myf.write('%6d %6d %6d \n'
                     % (self.kpt_n21, self.kpt_n22, self.kpt_n23))
            myf.write('%6d %6d %6d \n' 
                     % (self.kpt_n31, self.kpt_n32, self.kpt_n33))
            myf.write('  %10.6f %10.6f %10.6f \n' 
                     % (self.kpt_s1, self.kpt_s2, self.kpt_s3))
            myf.write('  }\n')

        #Dispersion parameters
        if (self.include_dispersion):
            myf.write('Dispersion = SlaterKirkwood {\n')
            myf.write(' PolarRadiusCharge = HybridDependentPol {\n')
            myf.write('\n')
            myf.write('  C={\n')
            myf.write('    CovalentRadius [Angstrom] = 0.8\n')
            myf.write('    HybridPolarisations [Angstrom^3,Angstrom,] = {\n')
            myf.write('      1.382 1.382 1.382 \
1.064 1.064 1.064 3.8 3.8 3.8 3.8 3.8 3.8 2.5\n')
            myf.write('    }\n')
            myf.write('  }\n')
            myf.write('\n')
            myf.write('  N={\n')
            myf.write('    CovalentRadius [Angstrom] = 0.8\n')
            myf.write('    HybridPolarisations [Angstrom^3,Angstrom,] = {\n')
            myf.write('      1.030 1.030 1.090 1.090 1.090 1.090 \
3.8 3.8 3.8 3.8 3.8 3.8 2.82\n')
            myf.write('    }\n')
            myf.write('  }\n')
            myf.write('\n')
            myf.write('  O={\n')
            myf.write('    CovalentRadius [Angstrom] = 0.8\n')
            myf.write('    HybridPolarisations [Angstrom^3,Angstrom,] = {\n')
            myf.write('    # All polarisabilities and radii set the same\n')
            myf.write('      0.560 0.560 0.000 0.000 0.000 0.000 \
3.8 3.8 3.8 3.8 3.8 3.8 3.15\n')
            myf.write('    }\n')
            myf.write('  }\n')
            myf.write('\n')


            myf.write('  H={\n')
            myf.write('    CovalentRadius [Angstrom] = 0.4\n')
            myf.write('    HybridPolarisations [Angstrom^3,Angstrom,] = {\n')
            myf.write('    # Different polarisabilities depending on the hybridisation\n')
            myf.write('      0.386 0.396 0.000 0.000 0.000 0.000 \
3.5 3.5 3.5 3.5 3.5 3.5 0.8\n')
            myf.write('    }\n')
            myf.write('  }\n')
            myf.write('\n')

            myf.write('  P={\n')
            myf.write('    CovalentRadius [Angstrom] = 0.9\n')
            myf.write('    HybridPolarisations [Angstrom^3,Angstrom,] = {\n')
            myf.write('    # Different polarisabilities depending \
on the hybridisation\n')
            myf.write('      1.600 1.600 1.600 1.600 1.600 1.600 \
4.7 4.7 4.7 4.7 4.7 4.7 4.50\n')
            myf.write('    }\n')
            myf.write('  }\n')
            myf.write('\n')

            myf.write('  S={\n')
            myf.write('    CovalentRadius [Angstrom] = 0.9\n')
            myf.write('    HybridPolarisations [Angstrom^3,Angstrom,] = {\n')
            myf.write('    # Different polarisabilities depending \
on the hybridisation\n')
            myf.write('      3.000 3.000 3.000 3.000 3.000 3.000 \
4.7 4.7 4.7 4.7 4.7 4.7 4.80\n')
            myf.write('    }\n')
            myf.write('  }\n')
            myf.write(' }\n')        
            myf.write('}\n') 
        myf.write('}\n') 
        myf.write('Options = {}\n')
        myf.write('ParserOptions = {\n')
        myf.write('  ParserVersion = 3\n')
        myf.write('}\n')
 
        myf.close()

    def change_atom_positions_dftb(self, atoms):
        """Write coordinates in DFTB+ input file dftb_in.hsd
        """

        filename = 'dftb_in.hsd'
        myf = open(filename)
        lines = myf.readlines()
        myf.close()

        myf = open(filename, 'w')

        coord = atoms.get_positions().copy()

        start_writing_coords = False
        stop_writing_coords = False
        i = 0
        for line in lines:
            if ('TypesAndCoordinates' in line):
                start_writing_coords = True
            if (start_writing_coords and not(stop_writing_coords)):
                if ('}' in line):
                    stop_writing_coords = True
            if (start_writing_coords and not(stop_writing_coords)and 
                not ('TypesAndCoordinates' in line)):
                atom_type_index = line.split()[0]
                myf.write('%6s  %20.14f  %20.14f  %20.14f \n'
                        % (atom_type_index,coord[i][0],coord[i][1],coord[i][2]))
                i = i + 1
            else:
                myf.write(line)

        myf.close()
