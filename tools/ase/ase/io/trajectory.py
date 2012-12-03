import os
import sys
import cPickle as pickle
import warnings
import errno

# pass for WindowsError on non-Win platforms
try:
    WindowsError
except NameError:
    class WindowsError(OSError): pass 

from ase.calculators.singlepoint import SinglePointCalculator
from ase.atoms import Atoms
from ase.parallel import rank, barrier
from ase.utils import devnull


class PickleTrajectory:
    """Reads/writes Atoms objects into a .traj file."""
    # Per default, write these quantities
    write_energy = True
    write_forces = True
    write_stress = True
    write_magmoms = True
    write_momenta = True
    write_info = True

    def __init__(self, filename, mode='r', atoms=None, master=None,
                 backup=True):
        """A PickleTrajectory can be created in read, write or append mode.

        Parameters:

        filename:
            The name of the parameter file.  Should end in .traj.

        mode='r':
            The mode.

            'r' is read mode, the file should already exist, and
            no atoms argument should be specified.

            'w' is write mode.  If the file already exists, is it
            renamed by appending .bak to the file name.  The atoms
            argument specifies the Atoms object to be written to the
            file, if not given it must instead be given as an argument
            to the write() method.

            'a' is append mode.  It acts a write mode, except that
            data is appended to a preexisting file.

        atoms=None:
            The Atoms object to be written in write or append mode.

        master=None:
            Controls which process does the actual writing. The
            default is that process number 0 does this.  If this
            argument is given, processes where it is True will write.

        backup=True:
            Use backup=False to disable renaming of an existing file.
        """

        self.numbers = None
        self.pbc = None
        self.sanitycheck = True
        self.pre_observers = []   # Callback functions before write
        self.post_observers = []  # Callback functions after write
        self.write_counter = 0    # Counter used to determine when callbacks
                                  # are called

        self.offsets = []
        if master is None:
            master = (rank == 0)
        self.master = master
        self.backup = backup
        self.set_atoms(atoms)
        self.open(filename, mode)

    def open(self, filename, mode):
        """Opens the file.

        For internal use only.
        """
        self.fd = filename
        if mode == 'r':
            if isinstance(filename, str):
                self.fd = open(filename, 'rb')
            self.read_header()
        elif mode == 'a':
            exists = True
            if isinstance(filename, str):
                exists = os.path.isfile(filename)
                if exists:
                    self.fd = open(filename, 'rb')
                    self.read_header()
                    self.fd.close()
                barrier()
                if self.master:
                    self.fd = open(filename, 'ab+')
                else:
                    self.fd = devnull
        elif mode == 'w':
            if self.master:
                if isinstance(filename, str):
                    if self.backup and os.path.isfile(filename):
                        try:
                            os.rename(filename, filename + '.bak')
                        except WindowsError, e:
                            # this must run on Win only! Not atomic!
                            if e.errno != errno.EEXIST:
                                raise
                            os.unlink(filename + '.bak')
                            os.rename(filename, filename + '.bak')
                    self.fd = open(filename, 'wb')
            else:
                self.fd = devnull
        else:
            raise ValueError('mode must be "r", "w" or "a".')

    def set_atoms(self, atoms=None):
        """Associate an Atoms object with the trajectory.

        Mostly for internal use.
        """
        if atoms is not None and not hasattr(atoms, 'get_positions'):
            raise TypeError('"atoms" argument is not an Atoms object.')
        self.atoms = atoms

    def read_header(self):
        self.fd.seek(0)
        try:
            if self.fd.read(len('PickleTrajectory')) != 'PickleTrajectory':
                raise IOError('This is not a trajectory file!')
            d = pickle.load(self.fd)
        except EOFError:
            raise EOFError('Bad trajectory file.')

        self.pbc = d['pbc']
        self.numbers = d['numbers']
        self.tags = d.get('tags')
        self.masses = d.get('masses')
        self.constraints = dict2constraints(d)
        self.offsets.append(self.fd.tell())

    def write(self, atoms=None):
        """Write the atoms to the file.

        If the atoms argument is not given, the atoms object specified
        when creating the trajectory object is used.
        """
        self._call_observers(self.pre_observers)
        if atoms is None:
            atoms = self.atoms

        if hasattr(atoms, 'interpolate'):
            # seems to be a NEB
            neb = atoms
            assert not neb.parallel
            try:
                neb.get_energies_and_forces(all=True)
            except AttributeError:
                pass
            for image in neb.images:
                self.write(image)
            return
        
        if len(self.offsets) == 0:
            self.write_header(atoms)
        else:
            if (atoms.pbc != self.pbc).any():
                raise ValueError('Bad periodic boundary conditions!')
            elif self.sanitycheck and len(atoms) != len(self.numbers):
                raise ValueError('Bad number of atoms!')
            elif self.sanitycheck and (atoms.numbers != self.numbers).any():
                raise ValueError('Bad atomic numbers!')
            
        if atoms.has('momenta'):
            momenta = atoms.get_momenta()
        else:
            momenta = None

        d = {'positions': atoms.get_positions(),
             'cell': atoms.get_cell(),
             'momenta': momenta}

        if atoms.get_calculator() is not None:
            if self.write_energy:
                d['energy'] = atoms.get_potential_energy()
            if self.write_forces:
                assert self.write_energy
                try:
                    d['forces'] = atoms.get_forces(apply_constraint=False)
                except NotImplementedError:
                    pass
            if self.write_stress:
                assert self.write_energy
                try:
                    d['stress'] = atoms.get_stress()
                except NotImplementedError:
                    pass

            if self.write_magmoms:
                try:
                    if atoms.calc.get_spin_polarized():
                        d['magmoms'] = atoms.get_magnetic_moments()
                except (NotImplementedError, AttributeError):
                    pass

        if 'magmoms' not in d and atoms.has('magmoms'):
            d['magmoms'] = atoms.get_initial_magnetic_moments()

        if self.write_info:
            d['info'] = stringnify_info(atoms.info)
            
        if self.master:
            pickle.dump(d, self.fd, protocol=-1)
        self.fd.flush()
        self.offsets.append(self.fd.tell())
        self._call_observers(self.post_observers)
        self.write_counter += 1

    def write_header(self, atoms):
        self.fd.write('PickleTrajectory')
        if atoms.has('tags'):
            tags = atoms.get_tags()
        else:
            tags = None
        if atoms.has('masses'):
            masses = atoms.get_masses()
        else:
            masses = None
        d = {'version': 3,
             'pbc': atoms.get_pbc(),
             'numbers': atoms.get_atomic_numbers(),
             'tags': tags,
             'masses': masses,
             'constraints': [],  # backwards compatibility
             'constraints_string': pickle.dumps(atoms.constraints)}
        pickle.dump(d, self.fd, protocol=-1)
        self.header_written = True
        self.offsets.append(self.fd.tell())

        # Atomic numbers and periodic boundary conditions are only
        # written once - in the header.  Store them here so that we can
        # check that they are the same for all images:
        self.numbers = atoms.get_atomic_numbers()
        self.pbc = atoms.get_pbc()
        
    def close(self):
        """Close the trajectory file."""
        self.fd.close()

    def __getitem__(self, i=-1):
        if isinstance(i, slice):
            return [self[j] for j in range(*i.indices(len(self)))]

        N = len(self.offsets)
        if 0 <= i < N:
            self.fd.seek(self.offsets[i])
            try:
                d = pickle.load(self.fd)
            except EOFError:
                raise IndexError
            if i == N - 1:
                self.offsets.append(self.fd.tell())
            try:
                magmoms = d['magmoms']
            except KeyError:
                magmoms = None
            atoms = Atoms(positions=d['positions'],
                          numbers=self.numbers,
                          cell=d['cell'],
                          momenta=d['momenta'],
                          magmoms=magmoms,
                          tags=self.tags,
                          masses=self.masses,
                          pbc=self.pbc,
                          info=unstringnify_info(d.get('info', {})),
                          constraint=[c.copy() for c in self.constraints])
            if 'energy' in d:
                calc = SinglePointCalculator(
                    d.get('energy', None), d.get('forces', None),
                    d.get('stress', None), magmoms, atoms)
                atoms.set_calculator(calc)
            return atoms

        if i >= N:
            for j in range(N - 1, i + 1):
                atoms = self[j]
            return atoms

        i = len(self) + i
        if i < 0:
            raise IndexError('Trajectory index out of range.')
        return self[i]

    def __len__(self):
        if len(self.offsets) == 0:
            return 0
        N = len(self.offsets) - 1
        while True:
            self.fd.seek(self.offsets[N])
            try:
                pickle.load(self.fd)
            except EOFError:
                return N
            self.offsets.append(self.fd.tell())
            N += 1

    def __iter__(self):
        del self.offsets[1:]
        return self

    def next(self):
        try:
            return self[len(self.offsets) - 1]
        except IndexError:
            raise StopIteration

    def guess_offsets(self):
        size = os.path.getsize(self.fd.name)

        while True:
            self.fd.seek(self.offsets[-1])
            try:
                pickle.load(self.fd)
            except:
                raise EOFError('Damaged trajectory file.')
            else:
                self.offsets.append(self.fd.tell())

            if self.offsets[-1] >= size:
                break

            if len(self.offsets) > 2:
                step1 = self.offsets[-1] - self.offsets[-2]
                step2 = self.offsets[-2] - self.offsets[-3]

                if step1 == step2:
                    m = int((size - self.offsets[-1]) / step1) - 1

                    while m > 1:
                        self.fd.seek(self.offsets[-1] + m * step1)
                        try:
                            pickle.load(self.fd)
                        except:
                            m = m / 2
                        else:
                            for i in range(m):
                                self.offsets.append(self.offsets[-1] + step1)
                            m = 0

    def pre_write_attach(self, function, interval=1, *args, **kwargs):
        """Attach a function to be called before writing begins.

        function: The function or callable object to be called.

        interval: How often the function is called.  Default: every time (1).

        All other arguments are stored, and passed to the function.
        """
        if not callable(function):
            raise ValueError('Callback object must be callable.')
        self.pre_observers.append((function, interval, args, kwargs))

    def post_write_attach(self, function, interval=1, *args, **kwargs):
        """Attach a function to be called after writing ends.

        function: The function or callable object to be called.

        interval: How often the function is called.  Default: every time (1).

        All other arguments are stored, and passed to the function.
        """
        if not callable(function):
            raise ValueError('Callback object must be callable.')
        self.post_observers.append((function, interval, args, kwargs))

    def _call_observers(self, obs):
        """Call pre/post write observers."""
        for function, interval, args, kwargs in obs:
            if self.write_counter % interval == 0:
                function(*args, **kwargs)


def stringnify_info(info):
    """Return a stringnified version of the dict *info* that is
    ensured to be picklable.  Items with non-string keys or
    unpicklable values are dropped and a warning is issued."""
    stringnified = {}
    for k, v in info.items():
        if not isinstance(k, basestring):
            warnings.warn('Non-string info-dict key is not stored in ' +
                          'trajectory: ' + repr(k), UserWarning)
            continue
        try:
            # Should highest protocol be used here for efficiency?
            # Protocol 2 seems not to raise an exception when one
            # tries to pickle a file object, so by using that, we
            # might end up with file objects in inconsistent states.
            s = pickle.dumps(v)
        except:
            warnings.warn('Skipping not picklable info-dict item: ' +
                          '"%s" (%s)' % (k, sys.exc_info()[1]), UserWarning)
        else:
            stringnified[k] = s
    return stringnified


def unstringnify_info(stringnified):
    """Convert the dict *stringnified* to a dict with unstringnified
    objects and return it.  Objects that cannot be unpickled will be
    skipped and a warning will be issued."""
    info = {}
    for k, s in stringnified.items():
        try:
            v = pickle.loads(s)
        except:
            warnings.warn('Skipping not unpicklable info-dict item: ' +
                          '"%s" (%s)' % (k, sys.exc_info()[1]), UserWarning)
        else:
            info[k] = v
    return info


def read_trajectory(filename, index=-1):
    traj = PickleTrajectory(filename, mode='r')

    if isinstance(index, int):
        return traj[index]
    else:
        # Here, we try to read only the configurations we need to read
        # and len(traj) should only be called if we need to as it will
        # read all configurations!

        # XXX there must be a simpler way?
        step = index.step or 1
        if step > 0:
            start = index.start or 0
            if start < 0:
                start += len(traj)
            stop = index.stop or len(traj)
            if stop < 0:
                stop += len(traj)
        else:
            if index.start is None:
                start = len(traj) - 1
            else:
                start = index.start
                if start < 0:
                    start += len(traj)
            if index.stop is None:
                stop = -1
            else:
                stop = index.stop
                if stop < 0:
                    stop += len(traj)
                    
        return [traj[i] for i in range(start, stop, step)]


def write_trajectory(filename, images):
    """Write image(s) to trajectory.

    Write also energy, forces, and stress if they are already
    calculated."""

    traj = PickleTrajectory(filename, mode='w')

    if not isinstance(images, (list, tuple)):
        images = [images]
        
    for atoms in images:
        # Avoid potentially expensive calculations:
        calc = atoms.get_calculator()
        if calc is not None:
            if  hasattr(calc, 'calculation_required'):
                if calc.calculation_required(atoms, ['energy']):
                    traj.write_energy = False
                if calc.calculation_required(atoms, ['forces']):
                    traj.write_forces = False
                if calc.calculation_required(atoms, ['stress']):
                    traj.write_stress = False
                if calc.calculation_required(atoms, ['magmoms']):
                    traj.write_magmoms = False
        else:
            traj.write_energy = False
            traj.write_forces = False
            traj.write_stress = False
            traj.write_magmoms = False
            
        traj.write(atoms)
    traj.close()


def dict2constraints(d):
    """Convert dict unpickled from trajectory file to list of constraints."""

    version = d.get('version', 1)

    if version == 1:
        return d['constraints']
    elif version in (2, 3):
        try:
            return pickle.loads(d['constraints_string'])
        except (AttributeError, KeyError, EOFError):
            warnings.warn('Could not unpickle constraints!')
            return []
    else:
        return []


def print_trajectory_info(filename):
    """Prints information about a PickleTrajectory file.

    Mainly intended to be called from a command line tool.
    """
    f = open(filename)
    hdr = 'PickleTrajectory'
    x = f.read(len(hdr))
    if x != hdr:
        raise ValueError('Not a PickleTrajectory file!')
    # Head header
    header = pickle.load(f)
    print('Header information of trajectory file %r:' % filename)
    print('  Version: %d' % header.get('version', 1))
    print('  Boundary conditions: %s' % header['pbc'])
    print('  Atomic numbers: shape = %s, type = %s' %
          (header['numbers'].shape, header['numbers'].dtype))
    if header.get('tags') is None:
        print('  Tags are absent.')
    else:
        print('  Tags: shape = %s, type = %s' %
              (header['tags'].shape, header['tags'].dtype))
    if header.get('masses') is None:
        print('  Masses are absent.')
    else:
        print('  Masses: shape = %s, type = %s' %
              (header['masses'].shape, header['masses'].dtype))
    constraints = dict2constraints(header)
    if constraints:
        print('  %d constraints are present.' % len(constraints))
    else:
        print('  No constraints.')

    after_header = f.tell()

    # Read the first frame
    frame = pickle.load(f)
    print('Contents of first frame:')
    for k, v in frame.items():
        if hasattr(v, 'shape'):
            print('  %s: shape = %s, type = %s' % (k, v.shape, v.dtype))
        else:
            print('  %s: %s' % (k, v))
    after_frame = f.tell()
    kB = 1024
    MB = 1024 * kB
    GB = 1024 * MB
    framesize = after_frame - after_header
    if framesize >= GB:
        print('Frame size: %.2f GB' % (1.0 * framesize / GB))
    elif framesize >= MB:
        print('Frame size: %.2f MB' % (1.0 * framesize / MB))
    else:
        print('Frame size: %.2f kB' % (1.0 * framesize / kB))

    # Print information about file size
    try:
        filesize = os.path.getsize(filename)
    except IOError:
        print('No information about the file size.')
    else:
        if filesize >= GB:
            print('File size: %.2f GB' % (1.0 * filesize / GB))
        elif filesize >= MB:
            print('File size: %.2f MB' % (1.0 * filesize / MB))
        else:
            print('File size: %.2f kB' % (1.0 * filesize / kB))
        
        nframes = (filesize - after_header) // framesize
        offset = nframes * framesize + after_header - filesize
        if offset == 0:
            if nframes == 1:
                print('Trajectory contains 1 frame.')
            else:
                print('Trajectory contains %d frames.' % nframes)
        else:
            print('Trajectory appears to contain approximately %d frames,' %
                  nframes)
            print('but the file size differs by %d bytes from the expected' %
                  -offset)
            print('value.')
