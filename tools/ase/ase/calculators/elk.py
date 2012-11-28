import os

import numpy as np

from ase.units import Bohr, Hartree

elk_parameters = {
    'swidth': Hartree,
    }

class ELK:
    def __init__(self, dir='.', xc=None, kpts=None, tasks=[0], **kwargs):
        for key, value in kwargs.items():
            if key in elk_parameters:
                kwargs[key] /= elk_parameters[key]
        try:
            self.exe = os.environ['ELK']
        except KeyError:
            self.exe = 'elk'

        if xc is not None:
            if 'xctype' in kwargs:
                raise ValueError("You can't use both 'xc' and 'xctype'!")
            else:
                kwargs['xctype'] = {'LDA': 3, # PW92
                                    'PBE': 20,
                                    'REVPBE': 21,
                                    'PBESOL': 22,
                                    'WC06': 26,
                                    'AM05': 30}[xc.upper()]

        if kpts is not None:
            if 'autokpt' in kwargs:
                if kwargs['autokpt']:
                    raise ValueError("You can't use both 'kpts' and 'autokpt'!")
            if 'ngridk' in kwargs:
                raise ValueError("You can't use both 'kpts' and 'ngridk'!")
            if 'vkloff' in kwargs:
                raise ValueError("You can't use both 'kpts' and 'vkloff'!")
            else:
                kwargs['ngridk'] = kpts
                vkloff = []
                for nk in kpts:
                    if nk % 2 == 0:  # shift kpoint away from gamma point
                        vkloff.append(0.5 / nk)
                    else:
                        vkloff.append(0)
                kwargs['vkloff'] = vkloff

        kwargs['tasks'] = tasks

        if 'rmt' in kwargs:
            self.rmt = kwargs['rmt']
            assert len(self.rmt.keys()) == len(list(set(self.rmt.keys()))), 'redundant rmt definitions'
        else:
            self.rmt = None

        self.parameters = kwargs
        if 'rmt' in self.parameters:
            self.parameters.pop('rmt') # this is not an elk keyword!

        self.dir = dir
        self.energy = None

        self.converged = False

    def update(self, atoms):
        if (not self.converged or
            len(self.numbers) != len(atoms) or
            (self.numbers != atoms.get_atomic_numbers()).any()):
            self.initialize(atoms)
            self.calculate(atoms)
        elif ((self.positions != atoms.get_positions()).any() or
              (self.pbc != atoms.get_pbc()).any() or
              (self.cell != atoms.get_cell()).any()):
            self.calculate(atoms)

    def initialize(self, atoms):
        self.numbers = atoms.get_atomic_numbers().copy()

        if not hasattr(self, 'spinpol'):
            self.spinpol = atoms.get_initial_magnetic_moments().any()

        self.write(atoms)

    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()

    def get_stress(self, atoms):
        self.update(atoms)
        return self.stress.copy()

    def get_ibz_k_points(self):
        return self.read_kpts(mode='ibz_k_points')

    def read_kpts(self, mode='ibz_k_points'):
        """ Returns list of kpts weights or kpts coordinates.  """
        values = []
        assert mode in ['ibz_k_points' , 'k_point_weights'], 'mode not in [\'ibz_k_points\' , \'k_point_weights\']'
        KPOINTS_file = '%s/KPOINTS.OUT' % self.dir
        if os.path.isfile(KPOINTS_file) or os.path.islink(KPOINTS_file):
            lines = open(KPOINTS_file).readlines()
            kpts = None
            for line in lines:
                if line.rfind(': nkpt') > -1:
                    kpts = int(line.split(':')[0].strip())
                    break
            assert not kpts is None
            text = lines[1:] # remove first line
            values = []
            for line in text:
                if mode == 'ibz_k_points':
                    b = [float(c.strip()) for c in line.split()[1:-3]]
                else:
                    b = [float(c.strip()) for c in line.split()[-2]]
                values.append(b)
            if len(values) == 0:
                values = None
            return np.array(values)

    def get_spin_polarized(self):
        return self.spinpol

    def get_number_of_spins(self):
        return 1 + int(self.spinpol)

    def get_magnetic_moment(self, atoms):
        self.update(atoms)
        return self.magnetic_moment

    def read_magnetic_moment(self):
        magmom = None
        INFO_file = '%s/INFO.OUT' % self.dir
        if os.path.isfile(INFO_file) or os.path.islink(INFO_file):
            lines = open(INFO_file).readlines()
            for line in lines:
                if line.rfind('total moment                :') > -1:
                    magmom = float(line.split(':')[1].strip()) # last iter
        return magmom

    def get_magnetic_moments(self, atoms):
        # not implemented yet, so
        # so set the total magnetic moment on the atom no. 0 and fill with 0.0
        self.update(atoms)
        magmoms = [0.0 for a in range(len(atoms))]
        magmoms[0] = self.get_magnetic_moment(atoms)
        return np.array(magmoms)

    def get_fermi_level(self):
        return self.read_fermi()

    def get_number_of_bands(self):
        return self.nbands

    def read_number_of_bands(self):
        nbands = None
        EIGVAL_file = '%s/EIGVAL.OUT' % self.dir
        if os.path.isfile(EIGVAL_file) or os.path.islink(EIGVAL_file):
            lines = open(EIGVAL_file).readlines()
            for line in lines:
                if line.rfind(': nstsv') > -1:
                    nbands = int(line.split(':')[0].strip())
                    break
        if self.get_spin_polarized():
            nbands = nbands / 2
        return nbands

    def get_number_of_iterations(self):
        return self.niter

    def read_number_of_iterations(self):
        niter = None
        INFO_file = '%s/INFO.OUT' % self.dir
        if os.path.isfile(INFO_file) or os.path.islink(INFO_file):
            lines = open(INFO_file).readlines()
            for line in lines:
                if line.rfind(' Loop number : ') > -1:
                    niter = int(line.split(':')[1].split()[0].strip()) # last iter
        return niter

    def get_electronic_temperature(self):
        return self.swidth

    def read_electronic_temperature(self):
        swidth = None
        INFO_file = '%s/INFO.OUT' % self.dir
        if os.path.isfile(INFO_file) or os.path.islink(INFO_file):
            text = open(INFO_file).read().lower()
            assert 'convergence targets achieved' in text
            assert not 'reached self-consistent loops maximum' in text
            for line in iter(text.split('\n')):
                if line.rfind('smearing width :') > -1:
                    swidth = float(line.split(':')[1].strip())
                    break
        return Hartree*swidth

    def get_number_of_electrons(self):
        return self.nelect

    def read_number_of_electrons(self):
        nelec = None
        INFO_file = '%s/INFO.OUT' % self.dir
        if os.path.isfile(INFO_file) or os.path.islink(INFO_file):
            text = open(INFO_file).read().lower()
            assert 'convergence targets achieved' in text
            assert not 'reached self-consistent loops maximum' in text
            # Total electronic charge
            for line in iter(text.split('\n')):
                if line.rfind('total electronic charge :') > -1:
                    nelec = float(line.split(':')[1].strip())
                    break
        return nelec

    def get_eigenvalues(self, kpt=0, spin=0):
        return self.read_eigenvalues(kpt, spin, 'eigenvalues')

    def get_occupation_numbers(self, kpt=0, spin=0):
        return self.read_eigenvalues(kpt, spin, 'occupations')

    def get_fermi_level(self):
        return self.fermi

    def calculate(self, atoms):
        self.positions = atoms.get_positions().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()

        self.initialize(atoms)

        assert os.system('cd %s&& %s ' % (self.dir, self.exe)) == 0
        self.read()

        self.converged = True

    def write(self, atoms):
        if not os.path.isdir(self.dir):
            os.mkdir(self.dir)
        fd = open('%s/elk.in' % self.dir, 'w')
        for key, value in self.parameters.items():
            fd.write('%s\n' % key)
            if isinstance(value, bool):
                fd.write('.%s.\n\n' % ('false', 'true')[value])
            elif isinstance(value, (int, float)):
                fd.write('%s\n\n' % value)
            else:
                fd.write('%s\n\n' % ' '.join([str(x) for x in value]))

        fd.write('avec\n')
        for vec in atoms.cell:
            fd.write('%.14f %.14f %.14f\n' % tuple(vec / Bohr))
        fd.write('\n')

        species = {}
        symbols = []
        for a, symbol in enumerate(atoms.get_chemical_symbols()):
            if symbol in species:
                species[symbol].append(a)
            else:
                species[symbol] = [a]
                symbols.append(symbol)
        fd.write('atoms\n%d\n' % len(species))
        scaled = atoms.get_scaled_positions()
        for symbol in symbols:
            fd.write("'%s.in'\n" % symbol)
            fd.write('%d\n' % len(species[symbol]))
            for a in species[symbol]:
                fd.write('%.14f %.14f %.14f 0.0 0.0 0.0\n' % tuple(scaled[a]))

        customspecies = self.rmt
        if customspecies:
            # custom species definitions
            fd.write("\n")
            sfile = os.path.join(os.environ['ELK_SPECIES_PATH'], 'elk.in')
            assert os.path.exists(sfile)
            slines = open(sfile, 'r').readlines()
            # remove unused species
            for s in customspecies.keys():
                if s not in species.keys():
                    customspecies.pop(s)
            # add undefined species with defaults
            for s in species.keys():
                if s not in customspecies.keys():
                    # use default rmt for undefined species
                    customspecies.update({s: 0.0})
            # write custom species into elk.in
            skeys = list(set(customspecies.keys())) # unique
            skeys.sort()
            for s in skeys:
                found = False
                for n, line in enumerate(slines):
                    if line.find("'" + s + "'") > -1:
                        begline = n - 1
                for n, line in enumerate(slines[begline:]):
                    if not line.strip(): # first empty line
                        endline = n
                        found = True
                        break
                assert found
                fd.write("species\n")
                # set rmt on third line
                rmt = customspecies[s]
                assert isinstance(rmt, (float,int))
                if rmt <= 0.0: # relative
                    # split needed because H is defined with comments
                    newrmt = float(slines[begline + 3].split()[0].strip()) + rmt
                else:
                    newrmt = rmt
                slines[begline + 3] = '%6s\n' % str(newrmt)
                for l in slines[begline: begline + endline]:
                    fd.write('%s' % l)
                fd.write("\n")
        else:
            # use default species
            # if sppath is present in elk.in it overwrites species blocks!
            fd.write("sppath\n'%s'\n\n" % os.environ['ELK_SPECIES_PATH'])


    def read_fermi(self):
        """Method that reads Fermi energy in Hartree from the output file
        and returns it in eV"""
        E_f=None
        INFO_file = '%s/INFO.OUT' % self.dir
        if os.path.isfile(INFO_file) or os.path.islink(INFO_file):
            text = open(INFO_file).read().lower()
            assert 'convergence targets achieved' in text
            assert not 'reached self-consistent loops maximum' in text
            for line in iter(text.split('\n')):
                if line.rfind('fermi                       :') > -1:
                    E_f = float(line.split(':')[1].strip())
            E_f = E_f*Hartree
        return E_f

    def read_eigenvalues(self, kpt=0, spin=0, mode='eigenvalues'):
        """ Returns list of last eigenvalues, occupations
        for given kpt and spin.  """
        values = []
        assert mode in ['eigenvalues' , 'occupations'], 'mode not in [\'eigenvalues\' , \'occupations\']'
        EIGVAL_file = '%s/EIGVAL.OUT' % self.dir
        if os.path.isfile(EIGVAL_file) or os.path.islink(EIGVAL_file):
            lines = open(EIGVAL_file).readlines()
            nstsv = None
            for line in lines:
                if line.rfind(': nstsv') > -1:
                    nstsv = int(line.split(':')[0].strip())
                    break
            assert not nstsv is None
            kpts = None
            for line in lines:
                if line.rfind(': nkpt') > -1:
                    kpts = int(line.split(':')[0].strip())
                    break
            assert not kpts is None
            text = lines[3:] # remove first 3 lines
            # find the requested k-point
            beg = 2 + (nstsv + 4) * kpt
            end = beg + nstsv
            if self.get_spin_polarized():
                # elk prints spin-up and spin-down together
                if spin == 0:
                    beg = beg
                    end = beg + nstsv / 2
                else:
                    beg = beg - nstsv / 2 - 3
                    end = end
            values = []
            for line in text[beg:end]:
                b = [float(c.strip()) for c in line.split()[1:]]
                values.append(b)
            if mode == 'eigenvalues':
                values = [Hartree*v[0] for v in values]
            else:
                values = [v[1] for v in values]
            if len(values) == 0:
                values = None
            return np.array(values)

    def read(self):
        fd = open('%s/TOTENERGY.OUT' % self.dir, 'r')
        self.energy = float(fd.readlines()[-1]) * Hartree
        # Forces:
        INFO_file = '%s/INFO.OUT' % self.dir
        if os.path.isfile(INFO_file) or os.path.islink(INFO_file):
            text = open(INFO_file).read().lower()
            assert 'convergence targets achieved' in text
            assert not 'reached self-consistent loops maximum' in text
            lines = iter(text.split('\n'))
            forces = []
            atomnum = 0
            for line in lines:
                if line.rfind('total force') > -1:
                    forces.append(np.array([float(f) for f in line.split(':')[1].split()]))
                    atomnum =+ 1
            self.forces = np.array(forces)
        else:
            raise RuntimeError
        # Stress
        self.stress = np.empty((3, 3))
        self.nbands = self.read_number_of_bands()
        self.nelect = self.read_number_of_electrons()
        self.swidth = self.read_electronic_temperature()
        self.niter = self.read_number_of_iterations()
        self.magnetic_moment = self.read_magnetic_moment()
