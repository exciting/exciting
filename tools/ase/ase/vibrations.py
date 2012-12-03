# -*- coding: utf-8 -*-

"""Vibrational modes."""

import pickle
from math import sin, pi, sqrt
from os import remove
from os.path import isfile
import sys

import numpy as np

import ase.units as units
from ase.io.trajectory import PickleTrajectory
from ase.parallel import rank, paropen
from ase.utils import opencew


class Vibrations:
    """Class for calculating vibrational modes using finite difference.

    The vibrational modes are calculated from a finite difference
    approximation of the Hessian matrix.

    The *summary()*, *get_energies()* and *get_frequencies()* methods all take
    an optional *method* keyword.  Use method='Frederiksen' to use the method
    described in:

      T. Frederiksen, M. Paulsson, M. Brandbyge, A. P. Jauho:
      "Inelastic transport theory from first-principles: methodology and
      applications for nanoscale devices", Phys. Rev. B 75, 205413 (2007)

    atoms: Atoms object
        The atoms to work on.
    indices: list of int
        List of indices of atoms to vibrate.  Default behavior is
        to vibrate all atoms.
    name: str
        Name to use for files.
    delta: float
        Magnitude of displacements.
    nfree: int
        Number of displacements per atom and cartesian coordinate, 2 and 4 are
        supported. Default is 2 which will displace each atom +delta and
        -delta for each cartesian coordinate.

    Example:

    >>> from ase import Atoms
    >>> from ase.calculators import EMT
    >>> from ase.optimize import BFGS
    >>> from ase.vibrations import Vibrations
    >>> n2 = Atoms('N2', [(0, 0, 0), (0, 0, 1.1)],
    ...            calculator=EMT())
    >>> BFGS(n2).run(fmax=0.01)
    BFGS:   0  16:01:21        0.440339       3.2518
    BFGS:   1  16:01:21        0.271928       0.8211
    BFGS:   2  16:01:21        0.263278       0.1994
    BFGS:   3  16:01:21        0.262777       0.0088
    >>> vib = Vibrations(n2)
    >>> vib.run()
    Writing vib.eq.pckl
    Writing vib.0x-.pckl
    Writing vib.0x+.pckl
    Writing vib.0y-.pckl
    Writing vib.0y+.pckl
    Writing vib.0z-.pckl
    Writing vib.0z+.pckl
    Writing vib.1x-.pckl
    Writing vib.1x+.pckl
    Writing vib.1y-.pckl
    Writing vib.1y+.pckl
    Writing vib.1z-.pckl
    Writing vib.1z+.pckl
    >>> vib.summary()
    ---------------------
    #    meV     cm^-1
    ---------------------
    0    0.0       0.0
    1    0.0       0.0
    2    0.0       0.0
    3    2.5      20.4
    4    2.5      20.4
    5  152.6    1230.8
    ---------------------
    Zero-point energy: 0.079 eV
    >>> vib.write_mode(-1)  # write last mode to trajectory file

    """

    def __init__(self, atoms, indices=None, name='vib', delta=0.01, nfree=2):
        assert nfree in [2, 4]
        self.atoms = atoms
        if indices is None:
            indices = range(len(atoms))
        self.indices = np.asarray(indices)
        self.name = name
        self.delta = delta
        self.nfree = nfree
        self.H = None
        self.ir = None

    def run(self):
        """Run the vibration calculations.

        This will calculate the forces for 6 displacements per atom +/-x,
        +/-y, +/-z. Only those calculations that are not already done will be
        started. Be aware that an interrupted calculation may produce an empty
        file (ending with .pckl), which must be deleted before restarting the
        job. Otherwise the forces will not be calculated for that
        displacement.

        Note that the calculations for the different displacements can be done
        simultaneously by several independent processes. This feature relies
        on the existence of files and the subsequent creation of the file in
        case it is not found.
        """

        filename = self.name + '.eq.pckl'
        fd = opencew(filename)
        if fd is not None:
            self.calculate(filename, fd)

        p = self.atoms.positions.copy()
        for a in self.indices:
            for i in range(3):
                for sign in [-1, 1]:
                    for ndis in range(1, self.nfree // 2 + 1):
                        filename = ('%s.%d%s%s.pckl' %
                                    (self.name, a, 'xyz'[i],
                                     ndis * ' +-'[sign]))
                        fd = opencew(filename)
                        if fd is not None:
                            disp = ndis * sign * self.delta
                            self.atoms.positions[a, i] = p[a, i] + disp
                            self.calculate(filename, fd)
                            self.atoms.positions[a, i] = p[a, i]

    def calculate(self, filename, fd):
        forces = self.atoms.get_forces()
        if self.ir:
            dipole = self.calc.get_dipole_moment(self.atoms)
        if rank == 0:
            if self.ir:
                pickle.dump([forces, dipole], fd)
                sys.stdout.write(
                    'Writing %s, dipole moment = (%.6f %.6f %.6f)\n' %
                    (filename, dipole[0], dipole[1], dipole[2]))
            else:
                pickle.dump(forces, fd)
                sys.stdout.write('Writing %s\n' % filename)
            fd.close()
        sys.stdout.flush()

    def clean(self):
        if isfile(self.name + '.eq.pckl'):
            remove(self.name + '.eq.pckl')

        for a in self.indices:
            for i in 'xyz':
                for sign in '-+':
                    for ndis in range(1, self.nfree // 2 + 1):
                        name = '%s.%d%s%s.pckl' % (self.name, a, i,
                                                   ndis * sign)
                        if isfile(name):
                            remove(name)

    def read(self, method='standard', direction='central'):
        self.method = method.lower()
        self.direction = direction.lower()
        assert self.method in ['standard', 'frederiksen']
        assert self.direction in ['central', 'forward', 'backward']

        n = 3 * len(self.indices)
        H = np.empty((n, n))
        r = 0
        if direction != 'central':
            feq = pickle.load(open(self.name + '.eq.pckl'))
        for a in self.indices:
            for i in 'xyz':
                name = '%s.%d%s' % (self.name, a, i)
                fminus = pickle.load(open(name + '-.pckl'))
                fplus = pickle.load(open(name + '+.pckl'))
                if self.method == 'frederiksen':
                    fminus[a] -= fminus.sum(0)
                    fplus[a] -= fplus.sum(0)
                if self.nfree == 4:
                    fminusminus = pickle.load(open(name + '--.pckl'))
                    fplusplus = pickle.load(open(name + '++.pckl'))
                    if self.method == 'frederiksen':
                        fminusminus[a] -= fminusminus.sum(0)
                        fplusplus[a] -= fplusplus.sum(0)
                if self.direction == 'central':
                    if self.nfree == 2:
                        H[r] = .5 * (fminus - fplus)[self.indices].ravel()
                    else:
                        H[r] = H[r] = (-fminusminus +
                                       8 * fminus -
                                       8 * fplus +
                                       fplusplus)[self.indices].ravel() / 12.0
                elif self.direction == 'forward':
                    H[r] = (feq - fplus)[self.indices].ravel()
                else:
                    assert self.direction == 'backward'
                    H[r] = (fminus - feq)[self.indices].ravel()
                H[r] /= 2 * self.delta
                r += 1
        H += H.copy().T
        self.H = H
        m = self.atoms.get_masses()
        if 0 in [m[index] for index in self.indices]:
            raise RuntimeError('Zero mass encountered in one or more of '
                               'the vibrated atoms. Use Atoms.set_masses()'
                               ' to set all masses to non-zero values.')

        self.im = np.repeat(m[self.indices] ** -0.5, 3)
        omega2, modes = np.linalg.eigh(self.im[:, None] * H * self.im)
        self.modes = modes.T.copy()

        # Conversion factor:
        s = units._hbar * 1e10 / sqrt(units._e * units._amu)
        self.hnu = s * omega2.astype(complex) ** 0.5

    def get_energies(self, method='standard', direction='central'):
        """Get vibration energies in eV."""

        if (self.H is None or method.lower() != self.method or
            direction.lower() != self.direction):
            self.read(method, direction)
        return self.hnu

    def get_frequencies(self, method='standard', direction='central'):
        """Get vibration frequencies in cm^-1."""

        s = 0.01 * units._e / units._c / units._hplanck
        return s * self.get_energies(method, direction)

    def summary(self, method='standard', direction='central', freq=None,
                log=sys.stdout):
        """Print a summary of the vibrational frequencies.

        Parameters:

        method : string
            Can be 'standard'(default) or 'Frederiksen'.
        direction: string
            Direction for finite differences. Can be one of 'central'
            (default), 'forward', 'backward'.
        freq : numpy array
            Optional. Can be used to create a summary on a set of known
            frequencies.
        log : if specified, write output to a different location than
            stdout. Can be an object with a write() method or the name of a
            file to create.
        """

        if isinstance(log, str):
            log = paropen(log, 'a')
        write = log.write

        s = 0.01 * units._e / units._c / units._hplanck
        if freq != None:
            hnu = freq / s
        else:
            hnu = self.get_energies(method, direction)
        write('---------------------\n')
        write('  #    meV     cm^-1\n')
        write('---------------------\n')
        for n, e in enumerate(hnu):
            if e.imag != 0:
                c = 'i'
                e = e.imag
            else:
                c = ' '
                e = e.real
            write('%3d %6.1f%s  %7.1f%s\n' % (n, 1000 * e, c, s * e, c))
        write('---------------------\n')
        write('Zero-point energy: %.3f eV\n' %
              self.get_zero_point_energy(freq=freq))

    def get_zero_point_energy(self, freq=None):
        if freq is None:
            return 0.5 * self.hnu.real.sum()
        else:
            s = 0.01 * units._e / units._c / units._hplanck
            return 0.5 * freq.real.sum() / s

    def get_mode(self, n):
        mode = np.zeros((len(self.atoms), 3))
        mode[self.indices] = (self.modes[n] * self.im).reshape((-1, 3))
        return mode

    def write_mode(self, n, kT=units.kB * 300, nimages=30):
        """Write mode to trajectory file."""
        mode = self.get_mode(n) * sqrt(kT / abs(self.hnu[n]))
        p = self.atoms.positions.copy()
        n %= 3 * len(self.indices)
        traj = PickleTrajectory('%s.%d.traj' % (self.name, n), 'w')
        calc = self.atoms.get_calculator()
        self.atoms.set_calculator()
        for x in np.linspace(0, 2 * pi, nimages, endpoint=False):
            self.atoms.set_positions(p + sin(x) * mode)
            traj.write(self.atoms)
        self.atoms.set_positions(p)
        self.atoms.set_calculator(calc)
        traj.close()

    def write_jmol(self):
        """Writes file for viewing of the modes with jmol."""

        fd = open(self.name + '.xyz', 'w')
        symbols = self.atoms.get_chemical_symbols()
        f = self.get_frequencies()
        for n in range(3 * len(self.indices)):
            fd.write('%6d\n' % len(self.atoms))
            if f[n].imag != 0:
                c = 'i'
                f[n] = f[n].imag
            else:
                c = ' '
            fd.write('Mode #%d, f = %.1f%s cm^-1' % (n, f[n], c))
            if self.ir:
                fd.write(', I = %.4f (D/Å)^2 amu^-1.\n' % self.intensities[n])
            else:
                fd.write('.\n')
            mode = self.get_mode(n)
            for i, pos in enumerate(self.atoms.positions):
                fd.write('%2s %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n' %
                         (symbols[i], pos[0], pos[1], pos[2],
                          mode[i, 0], mode[i, 1], mode[i, 2]))
        fd.close()
