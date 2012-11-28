"""This module defines an ASE interface to FLAPW code FLEUR.

http://www.flapw.de
"""

import os

from subprocess import Popen, PIPE

import re

import numpy as np

from ase.units import Hartree, Bohr

class FLEUR:
    """Class for doing FLEUR calculations.

    In order to use fleur one has to define the following environment
    variables:

    FLEUR_INPGEN path to the input generator (inpgen.x) of fleur

    FLEUR path to the fleur executable. Note that fleur uses different
    executable for real and complex cases (systems with/without inversion
    symmetry), so FLEUR must point to the correct executable.

    The initialize_density step can be performed in parallel
    only if run on one compute node. FLEUR_SERIAL is used for this step.

    It is probable that user needs to tune manually the input file before
    the actual calculation, so in addition to the standard
    get_potential_energy function this class defines the following utility
    functions:

    write_inp
        generate the input file `inp`
    initialize_density
        creates the initial density after possible manual edits of `inp`
    calculate
        convergence the total energy. With fleur, one specifies always
        only the number of SCF-iterations so this function launches
        the executable several times and monitors the convergence.
    relax
        Uses fleur's internal algorithm for structure
        optimization. Requires that the proper optimization parameters
        (atoms to optimize etc.) are specified by hand in `inp`

    """
    def __init__(self, xc='LDA', kpts=None, nbands=None, convergence=None,
                 width=None, kmax=None, mixer=None, maxiter=None,
                 maxrelax=20, workdir=None, equivatoms=True):

        """Construct FLEUR-calculator object.

        Parameters
        ==========
        xc: str
            Exchange-correlation functional. Must be one of LDA, PBE,
            RPBE.
        kpts: list of three int
            Monkhost-Pack sampling.
        nbands: int
            Number of bands. (not used at the moment)
        convergence: dictionary
            Convergence parameters (currently only energy in eV)
            {'energy' : float}
        width: float
            Fermi-distribution width in eV.
        kmax: float
            Plane wave cutoff in a.u. If kmax is set then:
            gmax = 3.0 * kmax
            gmaxxc = int(2.5 * kmax * 10)/10. (from set_inp.f)
        mixer: dictionary
            Mixing parameters imix, alpha, spinf
            {'imix' : int, 'alpha' : float, 'spinf' : float}
        maxiter: int
            Maximum number of SCF iterations (name in the code: itmax)
        maxrelax: int
            Maximum number of relaxation steps
        workdir: str
            Working directory for the calculation
        equivatoms: bool
            If False: generate inequivalent atoms (default is True).
            Setting to False allows one for example to calculate spin-polarized dimers.
            Ssee http://www.flapw.de/pm/index.php?n=User-Documentation.InputFileForTheInputGenerator.
        """

        self.xc = xc
        self.kpts = kpts
        self.nbands = nbands
        self.width = width
        self.kmax = kmax
        self.itmax_step_default = 9 # SCF steps per run (default)
        self.itmax_step = 5 # SCF steps per run
        assert self.itmax_step_default <= 9
        assert self.itmax_step <= self.itmax_step_default
        self.itmax_default = 40
        if maxiter is None:
            self.itmax = self.itmax_default
        else:
            self.itmax = maxiter
        self.maxrelax = maxrelax
        self.mixer = mixer

        if convergence:
            self.convergence = convergence 
            self.convergence['energy'] /= Hartree 
        else:
            self.convergence = {'energy' : 0.0001}

        self.start_dir = None
        self.workdir = workdir
        if self.workdir:
            self.start_dir = os.getcwd()
            if not os.path.isdir(workdir):
                os.mkdir(workdir)
        else:
            self.workdir = '.'
            self.start_dir = '.'

        self.equivatoms = equivatoms

        self.converged = False

    def run_executable(self, mode='fleur', executable='FLEUR'):

        assert executable in ['FLEUR', 'FLEUR_SERIAL']

        executable_use = executable
        if executable == 'FLEUR_SERIAL' and not os.environ.get(executable, ''):
            executable_use = 'FLEUR' # use FLEUR if FLEUR_SERIAL not set
        try:
            code_exe = os.environ[executable_use]
        except KeyError:
            raise RuntimeError('Please set ' + executable_use)
        p = Popen(code_exe, shell=True, stdin=PIPE, stdout=PIPE,
                  stderr=PIPE)
        stat = p.wait()
        out = p.stdout.read()
        err = p.stderr.read()
        print mode, ': stat= ', stat, ' out= ', out, ' err=', err
        # special handling of exit status from density generation and regular fleur.x
        if mode in ['density']:
            if '!' in err:
                raise RuntimeError(executable_use + ' exited with a code %s' % err)
        else:
            if stat != 0:
                raise RuntimeError(executable_use + ' exited with a code %d' % stat)


    def update(self, atoms):
        """Update a FLEUR calculation."""

        if (not self.converged or
            len(self.numbers) != len(atoms) or
            (self.numbers != atoms.get_atomic_numbers()).any()):
            self.initialize(atoms)
            self.calculate(atoms)
        elif ((self.positions != atoms.get_positions()).any() or
              (self.pbc != atoms.get_pbc()).any() or
              (self.cell != atoms.get_cell()).any()):
            self.converged = False
            self.initialize(atoms)
            self.calculate(atoms)

    def initialize(self, atoms):
        """Create an input file inp and generate starting density."""

        self.converged = False
        self.initialize_inp(atoms)
        self.initialize_density(atoms)

    def initialize_inp(self, atoms):
        """Create a inp file"""
        os.chdir(self.workdir)

        self.numbers = atoms.get_atomic_numbers().copy()
        self.positions = atoms.get_positions().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()

        # create the input
        self.write_inp(atoms)

        os.chdir(self.start_dir)

    def initialize_density(self, atoms):
        """Creates a new starting density."""

        os.chdir(self.workdir)
        # remove possible conflicting files
        files2remove = ['cdn1', 'fl7para', 'stars', 'wkf2', 'enpara',
                        'kpts', 'broyd', 'broyd.7', 'tmat', 'tmas']
        if 0:
            # avoid STOP bzone3 error by keeping the kpts file
            files2remove.remove('kpts')

        for f in files2remove:
            if os.path.isfile(f):
                os.remove(f)

        # generate the starting density
        os.system("sed -i -e 's/strho=./strho=T/' inp")
        self.run_executable(mode='density', executable='FLEUR_SERIAL')
        os.system("sed -i -e 's/strho=./strho=F/' inp")

        os.chdir(self.start_dir)
        # generate spin-polarized density
        # http://www.flapw.de/pm/index.php?n=User-Documentation.Magnetism
        if atoms.get_initial_magnetic_moments().sum() > 0.0:
            os.chdir(self.workdir)
            # generate cdnc file (1 SCF step: swsp=F - non-magnetic)
            os.system("sed -i -e 's/itmax=.*,maxiter/itmax= 1,maxiter/' inp")
            self.run_executable(mode='cdnc', executable='FLEUR')
            sedline = "'s/itmax=.*,maxiter/itmax= '"
            sedline += str(self.itmax_step_default) + "',maxiter/'"
            os.system("sed -i -e " + sedline + " inp")
            # generate spin polarized density (swsp=T)
            os.system("sed -i -e 's/swsp=./swsp=T/' inp")
            self.run_executable(mode='swsp', executable='FLEUR_SERIAL')
            # restore swsp=F
            os.system("sed -i -e 's/swsp=./swsp=F/' inp")
            os.chdir(self.start_dir)

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)

        if force_consistent:
            return self.efree * Hartree
        else:
            # Energy extrapolated to zero Kelvin:
            return  (self.etotal + self.efree) / 2 * Hartree

    def get_number_of_iterations(self, atoms):
        self.update(atoms)
        return self.niter

    def get_forces(self, atoms):
        self.update(atoms)
        # electronic structure is converged, so let's calculate forces:
        # TODO
        return np.array((0.0, 0.0, 0.0))

    def get_stress(self, atoms):
        raise NotImplementedError

    def get_dipole_moment(self, atoms):
        """Returns total dipole moment of the system."""
        raise NotImplementedError

    def calculate(self, atoms):
        """Converge a FLEUR calculation to self-consistency.

           Input files should be generated before calling this function
           FLEUR performs always fixed number of SCF steps. This function
           reduces the number of iterations gradually, however, a minimum
           of five SCF steps is always performed.
        """

        os.chdir(self.workdir)

        self.niter = 0
        out = ''
        err = ''
        while not self.converged:
            if self.niter > self.itmax:
                raise RuntimeError('FLEUR failed to convergence in %d iterations' % self.itmax)

            self.run_executable(mode='fleur', executable='FLEUR')

            # catenate new output with the old one
            os.system('cat out >> out.old')
            self.read()
            self.check_convergence()

        if os.path.exists('out.old'): os.rename('out.old', 'out')
        # After convergence clean up broyd* files
        os.system('rm -f broyd*')
        os.chdir(self.start_dir)
        return out, err

    def relax(self, atoms):
        """Currently, user has to manually define relaxation parameters
           (atoms to relax, relaxation directions, etc.) in inp file
           before calling this function."""

        nrelax = 0
        relaxed = False
        while not relaxed:
            # Calculate electronic structure
            self.calculate(atoms)
            # Calculate the Pulay forces
            os.system("sed -i -e 's/l_f=./l_f=T/' inp")
            while True:
                self.converged = False
                out, err = self.calculate(atoms)
                if 'GEO new' in err:
                    os.chdir(self.workdir)
                    os.rename('inp_new', 'inp')
                    os.chdir(self.start_dir)
                    break
            if 'GEO: Des woas' in err:
                relaxed = True
                break
            nrelax += 1
            # save the out and cdn1 files
            os.system('cp out out_%d' % nrelax)
            os.system('cp cdn1 cdn1_%d' % nrelax)
            if nrelax > self.maxrelax:
                raise RuntimeError('Failed to relax in %d iterations' % self.maxrelax)
            self.converged = False


    def write_inp(self, atoms):
        """Write the `inp` input file of FLEUR.

        First, the information from Atoms is written to the simple input
        file and the actual input file `inp` is then generated with the
        FLEUR input generator. The location of input generator is specified
        in the environment variable FLEUR_INPGEN.

        Finally, the `inp` file is modified according to the arguments of
        the FLEUR calculator object.
        """

        fh = open('inp_simple', 'w')
        fh.write('FLEUR input generated with ASE\n')
        fh.write('\n')

        if atoms.pbc[2]:
            film = 'f'
        else:
            film = 't'
        fh.write('&input film=%s /' % film)
        fh.write('\n')

        for vec in atoms.get_cell():
            fh.write(' ')
            for el in vec:
                fh.write(' %21.16f' % (el/Bohr))
            fh.write('\n')
        fh.write(' %21.16f\n' % 1.0)
        fh.write(' %21.16f %21.16f %21.16f\n' % (1.0, 1.0, 1.0))
        fh.write('\n')

        natoms = len(atoms)
        fh.write(' %6d\n' % natoms)
        positions = atoms.get_scaled_positions()
        if not atoms.pbc[2]:
            # in film calculations z position has to be in absolute
            # coordinates and symmetrical
            cart_pos = atoms.get_positions()
            cart_pos[:, 2] -= atoms.get_cell()[2, 2]/2.0
            positions[:, 2] = cart_pos[:, 2] / Bohr
        atomic_numbers = atoms.get_atomic_numbers()
        for n, (Z, pos) in enumerate(zip(atomic_numbers, positions)):
            if self.equivatoms:
                fh.write('%3d' % Z)
            else:
                # generate inequivalent atoms, by using non-integer Z
                # (only the integer part will be used as Z of the atom)
                # see http://www.flapw.de/pm/index.php?n=User-Documentation.InputFileForTheInputGenerator
                fh.write('%3d.%04d' % (Z, n)) # MDTMP don't think one can calculate more that 10**4 atoms
            for el in pos:
                fh.write(' %21.16f' % el)
            fh.write('\n')

        # avoid "STOP read_record: ERROR reading input"
        fh.write('&end /')

        fh.close()
        try:
            inpgen = os.environ['FLEUR_INPGEN']
        except KeyError:
            raise RuntimeError('Please set FLEUR_INPGEN')

        # rename the previous inp if it exists
        if os.path.isfile('inp'):
            os.rename('inp', 'inp.bak')
        os.system('%s < inp_simple' % inpgen)

        # read the whole inp-file for possible modifications
        fh = open('inp', 'r')
        lines = fh.readlines()
        fh.close()


        window_ln = -1
        for ln, line in enumerate(lines):
            # XC potential
            if line.startswith('pbe'):
                if self.xc == 'PBE':
                    pass
                elif self.xc == 'RPBE':
                    lines[ln] = 'rpbe   non-relativi\n'
                elif self.xc == 'LDA':
                    lines[ln] = 'mjw    non-relativic\n'
                    del lines[ln+1]
                else:
                    raise RuntimeError('XC-functional %s is not supported' % self.xc)
            if line.startswith('Window'):
                # few things are set around this line
                window_ln = ln
            # kmax
            if self.kmax and ln == window_ln:
                line = '%10.5f\n' % self.kmax
                lines[ln+2] = line
            # gmax   cutoff for PW-expansion of potential & density  ( > 2*kmax)
            # gmaxxc cutoff for PW-expansion of XC-potential ( > 2*kmax, < gmax)
            if self.kmax and line.startswith('vchk'):
                gmax = 3. * self.kmax
                line = ' %10.6f %10.6f\n' % (gmax, int(2.5 * self.kmax * 10)/10.)
                lines[ln-1] = line
            # Fermi width
            if self.width and line.startswith('gauss'):
                line = 'gauss=F   %7.5ftria=F\n' % (self.width / Hartree)
                lines[ln] = line
            # kpts
            if self.kpts and line.startswith('nkpt'):
                line = 'nkpt=      nx=%2d,ny=%2d,nz=%2d\n' % (self.kpts[0],
                                                              self.kpts[1],
                                                              self.kpts[2])
                lines[ln] = line
            # itmax
            if self.itmax < self.itmax_step_default and line.startswith('itmax'):
                # decrease number of SCF steps; increasing is done by 'while not self.converged:'
                lsplit = line.split(',')
                if lsplit[0].find('itmax') != -1:
                    lsplit[0] = 'itmax=' + ('%2d' % self.itmax)
                    lines[ln] = ",".join(lsplit)
            # Mixing
            if self.mixer and line.startswith('itmax'):
                imix = self.mixer['imix']
                alpha = self.mixer['alpha']
                spinf = self.mixer['spinf']
                line_end = 'imix=%2d,alpha=%6.2f,spinf=%6.2f\n' % (imix,
                                                                   alpha,
                                                                   spinf)
                line = line[:21] + line_end
                lines[ln] = line
            # jspins and swsp
            if atoms.get_initial_magnetic_moments().sum() > 0.0:
                assert not self.equivatoms, 'equivatoms currently not allowed in magnetic systems'
                if line.find('jspins=1') != -1:
                    lines[ln] = line.replace('jspins=1', 'jspins=2')
                if line.startswith('swsp=F'):
                    # setting initial magnetic moments for all atom types
                    lines[ln] = 'swsp=F'
                    for m in atoms.get_initial_magnetic_moments():
                        lines[ln] += (' %5.2f' % m)
                    lines[ln] += '\n'
            # inpgen produces incorrect symbol 'J' for Iodine
            if line.startswith(' J  53'):
                lines[ln] = lines[ln].replace(' J  53', ' I  53')

        # write everything back to inp
        fh = open('inp', 'w')
        for line in lines:
            fh.write(line)
        fh.close()

    def read(self):
        """Read results from FLEUR's text-output file `out`."""

        lines = open('out', 'r').readlines()

        # total energies
        self.total_energies = []
        pat = re.compile('(.*total energy=)(\s)*([-0-9.]*)')
        for line in lines:
            m = pat.match(line)
            if m:
                self.total_energies.append(float(m.group(3)))
        self.etotal = self.total_energies[-1]

        # free_energies
        self.free_energies = []
        pat = re.compile('(.*free energy=)(\s)*([-0-9.]*)')
        for line in lines:
            m = pat.match(line)
            if m:
                self.free_energies.append(float(m.group(3)))
        self.efree = self.free_energies[-1]

        # TODO forces, charge density difference...

    def check_convergence(self):
        """Check the convergence of calculation"""
        energy_error = np.ptp(self.total_energies[-3:])
        self.converged = energy_error < self.convergence['energy']

        # TODO check charge convergence

        # reduce the itmax in inp
        lines = open('inp', 'r').readlines()
        pat = re.compile('(itmax=)([ 0-9]*)')
        fh = open('inp', 'w')
        for line in lines:
            m = pat.match(line)
            if m:
                itmax = int(m.group(2))
                self.niter += itmax
                itmax_new = itmax / 2
                itmax = max(self.itmax_step, itmax_new)
                line = 'itmax=%2d' % itmax + line[8:]
            fh.write(line)
        fh.close()
