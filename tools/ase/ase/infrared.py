
# -*- coding: utf-8 -*-

"""Infrared intensities"""

import pickle
from math import sin, pi, sqrt, exp, log

import numpy as np

import ase.units as units
from ase.io.trajectory import PickleTrajectory
from ase.parallel import rank, barrier, parprint
from ase.vibrations import Vibrations


class InfraRed(Vibrations):
    """Class for calculating vibrational modes and infrared intensities
    using finite difference.

    The vibrational modes are calculated from a finite difference
    approximation of the Dynamical matrix and the IR intensities from
    a finite difference approximation of the gradient of the dipole
    moment. The method is described in:

      D. Porezag, M. R. Pederson:
      "Infrared intensities and Raman-scattering activities within
      density-functional theory",
      Phys. Rev. B 54, 7830 (1996)

    The calculator object (calc) linked to the Atoms object (atoms) must 
    have the attribute:
    
    >>> calc.get_dipole_moment(atoms)

    In addition to the methods included in the ``Vibrations`` class
    the ``InfraRed`` class introduces two new methods;
    *get_spectrum()* and *write_spectra()*. The *summary()*, *get_energies()*, 
    *get_frequencies()*, *get_spectrum()* and *write_spectra()*
    methods all take an optional *method* keyword.  Use
    method='Frederiksen' to use the method described in:

      T. Frederiksen, M. Paulsson, M. Brandbyge, A. P. Jauho:
      "Inelastic transport theory from first-principles: methodology
      and applications for nanoscale devices", 
      Phys. Rev. B 75, 205413 (2007) 

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
        Number of displacements per degree of freedom, 2 or 4 are
        supported. Default is 2 which will displace each atom +delta
        and -delta in each cartesian direction.
    directions: list of int
        Cartesian coordinates to calculate the gradient of the dipole moment in. 
        For example directions = 2 only dipole moment in the z-direction will
        be considered, whereas for directions = [0, 1] only the dipole
        moment in the xy-plane will be considered. Default behavior is to
        use the dipole moment in all directions.

    Example:
    
    >>> from ase.io import read
    >>> from ase.calculators.vasp import Vasp
    >>> from ase.infrared import InfraRed
    >>> water = read('water.traj')  # read pre-relaxed structure of water molecule
    >>> calc = Vasp(prec='Accurate',
    ...             ediff=1E-8,
    ...             isym=0,
    ...             idipol=4,       # calculate the total dipole moment
    ...             dipol=water.get_center_of_mass(scaled=True),
    ...             ldipol=True)
    >>> water.set_calculator(calc)
    >>> ir = InfraRed(water)
    >>> ir.run()
    >>> ir.summary()
    -------------------------------------
    Mode    Frequency        Intensity
    #    meV     cm^-1   (D/Å)^2 amu^-1
    -------------------------------------
    0   16.9i    136.2i     1.6108
    1   10.5i     84.9i     2.1682
    2    5.1i     41.1i     1.7327
    3    0.3i      2.2i     0.0080
    4    2.4      19.0      0.1186
    5   15.3     123.5      1.4956
    6  195.5    1576.7      1.6437
    7  458.9    3701.3      0.0284
    8  473.0    3814.6      1.1812
    -------------------------------------
    Zero-point energy: 0.573 eV
    Static dipole moment: 1.833 D
    Maximum force on atom in `equilibrium`: 0.0026 eV/Å



    This interface now also works for calculator 'siesta', 
    (added get_dipole_moment for siesta).

    Example:

    >>> #!/usr/bin/env python

    >>> from ase.io import read
    >>> from ase.calculators.siesta import Siesta
    >>> from ase.infrared import InfraRed

    >>> bud = read('bud1.xyz')

    >>> calc = Siesta(label='bud',
    ...       meshcutoff=250 * Ry,
    ...       basis='DZP',
    ...       kpts=[1, 1, 1])

    >>> calc.set_fdf('DM.MixingWeight', 0.08)
    >>> calc.set_fdf('DM.NumberPulay', 3)
    >>> calc.set_fdf('DM.NumberKick', 20)
    >>> calc.set_fdf('DM.KickMixingWeight', 0.15)
    >>> calc.set_fdf('SolutionMethod',      'Diagon')
    >>> calc.set_fdf('MaxSCFIterations', 500)
    >>> calc.set_fdf('PAO.BasisType',  'split')
    >>> #50 meV = 0.003674931 * Ry
    >>> calc.set_fdf('PAO.EnergyShift', 0.003674931 * Ry )
    >>> calc.set_fdf('LatticeConstant', 1.000000 * Ang)
    >>> calc.set_fdf('WriteCoorXmol',       'T')

    >>> bud.set_calculator(calc)

    >>> ir = InfraRed(bud)
    >>> ir.run()
    >>> ir.summary()




    """
    def __init__(self, atoms, indices=None, name='ir', delta=0.01, nfree=2, directions=None):
        assert nfree in [2, 4]
        self.atoms = atoms
        if atoms.constraints:
            print "WARNING! \n Your Atoms object is constrained. Some forces may be unintended set to zero. \n"
        self.calc = atoms.get_calculator()
        if indices is None:
            indices = range(len(atoms))
        self.indices = np.asarray(indices)
        self.nfree = nfree
        self.name = name+'-d%.3f' % delta
        self.delta = delta
        self.H = None
        if directions is None:
            self.directions = np.asarray([0, 1, 2])
        else:
            self.directions = np.asarray(directions)
        self.ir = True

    def read(self, method='standard', direction='central'):
        self.method = method.lower()
        self.direction = direction.lower()
        assert self.method in ['standard', 'frederiksen']
        if direction != 'central':
            raise NotImplementedError('Only central difference is implemented at the moment.')

        # Get "static" dipole moment and forces
        name = '%s.eq.pckl' % self.name
        [forces_zero, dipole_zero] = pickle.load(open(name))
        self.dipole_zero = (sum(dipole_zero**2)**0.5) / units.Debye
        self.force_zero = max([sum((forces_zero[j])**2)**0.5 for j in self.indices])

        ndof = 3 * len(self.indices)
        H = np.empty((ndof, ndof))
        dpdx = np.empty((ndof, 3))
        r = 0
        for a in self.indices:
            for i in 'xyz':
                name = '%s.%d%s' % (self.name, a, i)
                [fminus, dminus] = pickle.load(open(name + '-.pckl'))
                [fplus, dplus] = pickle.load(open(name + '+.pckl'))
                if self.nfree == 4:
                    [fminusminus, dminusminus] = pickle.load(open(name + '--.pckl'))
                    [fplusplus, dplusplus] = pickle.load(open(name + '++.pckl'))
                if self.method == 'frederiksen':
                    fminus[a] += -fminus.sum(0)
                    fplus[a] += -fplus.sum(0)
                    if self.nfree == 4:
                        fminusminus[a] += -fminus.sum(0)
                        fplusplus[a] += -fplus.sum(0)
                if self.nfree == 2:
                    H[r] = (fminus - fplus)[self.indices].ravel() / 2.0
                    dpdx[r] = (dminus - dplus)
                if self.nfree == 4:
                    H[r] = (-fminusminus+8*fminus-8*fplus+fplusplus)[self.indices].ravel() / 12.0
                    dpdx[r] = (-dplusplus + 8*dplus - 8*dminus +dminusminus) / 6.0
                H[r] /= 2 * self.delta
                dpdx[r] /= 2 * self.delta
                for n in range(3):
                    if n not in self.directions:
                        dpdx[r][n] = 0
                        dpdx[r][n] = 0
                r += 1
        # Calculate eigenfrequencies and eigenvectors
        m = self.atoms.get_masses()
        H += H.copy().T
        self.H = H
        m = self.atoms.get_masses()
        self.im = np.repeat(m[self.indices]**-0.5, 3)
        omega2, modes = np.linalg.eigh(self.im[:, None] * H * self.im)
        self.modes = modes.T.copy()

        # Calculate intensities
        dpdq = np.array([dpdx[j]/sqrt(m[self.indices[j/3]]*units._amu/units._me) for j in range(ndof)])
        dpdQ = np.dot(dpdq.T, modes)
        dpdQ = dpdQ.T
        intensities = np.array([sum(dpdQ[j]**2) for j in range(ndof)])
        # Conversion factor:
        s = units._hbar * 1e10 / sqrt(units._e * units._amu)
        self.hnu = s * omega2.astype(complex)**0.5
        # Conversion factor from atomic units to (D/Angstrom)^2/amu.
        conv = (1.0 / units.Debye)**2*units._amu/units._me
        self.intensities = intensities*conv

    def summary(self, method='standard', direction='central', 
                intensity_unit='(D/A)2/amu'):
        hnu = self.get_energies(method, direction)
        s = 0.01 * units._e / units._c / units._hplanck
        if intensity_unit == '(D/A)2/amu':
            iu = 1.0
            iu_string = '(D/Å)^2 amu^-1'
            iu_format = '%9.4f'
        elif intensity_unit == 'km/mol':
            # conversion factor from Porezag PRB 54 (1996) 7830
            iu = 42.255
            iu_string = '   km/mol'
            iu_format = ' %7.1f'
        else:
            raise RuntimeError('Intensity unit >' + intensity_unit +
                               '< unknown.')
        parprint('-------------------------------------')
        parprint(' Mode    Frequency        Intensity')
        parprint('  #    meV     cm^-1   ' + iu_string)
        parprint('-------------------------------------')
        for n, e in enumerate(hnu):
            if e.imag != 0:
                c = 'i'
                e = e.imag
            else:
                c = ' '
            parprint(('%3d %6.1f%s  %7.1f%s  ' + iu_format) % 
                     (n, 1000 * e, c, s * e, c, iu * self.intensities[n]))
        parprint('-------------------------------------')
        parprint('Zero-point energy: %.3f eV' % self.get_zero_point_energy())
        parprint('Static dipole moment: %.3f D' % self.dipole_zero)
        parprint('Maximum force on atom in `equilibrium`: %.4f eV/Å' % 
                  self.force_zero)
        parprint()

    def get_spectrum(self, start=800, end=4000, npts=None, width=4, type='Gaussian', method='standard', direction='central'):
        """Get infrared spectrum.

        The method returns wavenumbers in cm^-1 with corresonding absolute infrared intensity.
        Start and end point, and width of the Gaussian/Lorentzian should be given in cm^-1."""

        self.type = type.lower()
        assert self.type in ['gaussian', 'lorentzian']
        if not npts: 
            npts = (end-start)/width*10+1
        frequencies = self.get_frequencies(method, direction).real
        intensities=self.intensities
        if type == 'lorentzian':
            intensities = intensities*width*pi/2.
        else:
            sigma = width/2./sqrt(2.*log(2.))
        #Make array with spectrum data
        spectrum = np.empty(npts,np.float)
        energies = np.empty(npts,np.float)
        ediff = (end-start)/float(npts-1)
        energies = np.arange(start, end+ediff/2, ediff)
        for i, energy in enumerate(energies):
            energies[i] = energy
            if type == 'lorentzian':
                spectrum[i] = (intensities*0.5*width/pi/((frequencies-energy)**2+0.25*width**2)).sum()
            else:
                spectrum[i] = (intensities*np.exp(-(frequencies - energy)**2/2./sigma**2)).sum()
        return [energies, spectrum]

    def write_spectra(self, out='ir-spectra.dat', start=800, end=4000, npts=None, width=10, type='Gaussian', method='standard', direction='central'):
        """Write out infrared spectrum to file.

        First column is the wavenumber in cm^-1, the second column the absolute infrared intensities, and
        the third column the absorbance scaled so that data runs from 1 to 0. Start and end 
        point, and width of the Gaussian/Lorentzian should be given in cm^-1."""
        energies, spectrum = self.get_spectrum(start, end, npts, width, type, method, direction)

        #Write out spectrum in file. First column is absolute intensities. 
        #Second column is absorbance scaled so that data runs from 1 to 0
        spectrum2 = 1. - spectrum/spectrum.max()
        outdata = np.empty([len(energies), 3])
        outdata.T[0] = energies
        outdata.T[1] = spectrum
        outdata.T[2] = spectrum2
        fd = open(out, 'w')
        for row in outdata:
            fd.write('%.3f  %15.5e  %15.5e \n' % (row[0], row[1], row[2]) )
        fd.close()
        #np.savetxt(out, outdata, fmt='%.3f  %15.5e  %15.5e')
