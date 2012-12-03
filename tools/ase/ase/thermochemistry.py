#!/usr/bin/env python
"""Modules for calculating thermochemical information from computational
outputs."""

import os
import sys
import numpy as np

from ase import units


def rotationalinertia(atoms):
    """Calculates the three principle moments of inertia for an ASE atoms
    object. This uses the atomic masses from ASE, which (if not explicitly
    specified by the user) gives an inexact approximation of an isotopically
    averaged result. Units are in amu*angstroms**2."""

    # Calculate the center of mass.
    xcm, ycm, zcm = atoms.get_center_of_mass()
    masses = atoms.get_masses()

    # Calculate moments of inertia in the current frame of reference.
    Ixx = 0.
    Iyy = 0.
    Izz = 0.
    Ixy = 0.
    Ixz = 0.
    Iyz = 0.
    for index, atom in enumerate(atoms):
        m = masses[index]
        x = atom.x - xcm
        y = atom.y - ycm
        z = atom.z - zcm
        Ixx += m * (y**2. + z**2.)
        Iyy += m * (x**2. + z**2.)
        Izz += m * (x**2. + y**2.)
        Ixy += m * x * y
        Ixz += m * x * z
        Iyz += m * y * z
    # Create the inertia tensor in the current frame of reference.
    I_ = np.matrix([[ Ixx, -Ixy, -Ixz],
                    [-Ixy,  Iyy, -Iyz],
                    [-Ixz, -Iyz,  Izz]])
    # Find the eigenvalues, which are the principle moments of inertia.
    I = np.linalg.eigvals(I_)
    return I


class ThermoChem:
    """Base class containing common methods used in thermochemistry
    calculations."""

    def get_ZPE_correction(self):
        """Returns the zero-point vibrational energy correction in eV."""
        zpe = 0.
        for energy in self.vib_energies:
            zpe += 0.5 * energy
        return zpe

    def _vibrational_energy_contribution(self, temperature):
        """Calculates the change in internal energy due to vibrations from
        0K to the specified temperature for a set of vibrations given in
        eV and a temperature given in Kelvin. Returns the energy change
        in eV."""
        kT = units.kB * temperature
        dU = 0.
        for energy in self.vib_energies:
            dU += energy / (np.exp(energy/kT) - 1.)
        return dU

    def _vibrational_entropy_contribution(self, temperature):
        """Calculates the entropy due to vibrations for a set of vibrations
        given in eV and a temperature given in Kelvin.  Returns the entropy
        in eV/K."""
        kT = units.kB * temperature
        S_v = 0.
        for energy in self.vib_energies:
            x = energy / kT
            S_v += x / (np.exp(x)-1.) - np.log(1-np.exp(-x))
        S_v *= units.kB
        return S_v

    def _vprint(self, text):
        """Print output if verbose flag True."""
        if self.verbose:
            sys.stdout.write(text + os.linesep)


class HarmonicThermo(ThermoChem):
    """Class for calculating thermodynamic properties in the approximation
    that all degrees of freedom are treated harmonically. Often used for
    adsorbates.
    
    Inputs:

    vib_energies : list
        a list of the harmonic energies of the adsorbate (e.g., from
        ase.vibrations.Vibrations.get_energies). The number of
        energies should match the number of degrees of freedom of the
        adsorbate; i.e., 3*n, where n is the number of atoms. Note that
        this class does not check that the user has supplied the correct
        number of energies. Units of energies are eV.
    electronicenergy : float
        the electronic energy in eV
        (If the electronicenergy is unspecified, then the methods of this
        class can be interpreted as the energy corrections.)
    """

    def __init__(self, vib_energies, electronicenergy=None):
        self.vib_energies = vib_energies
        # Check for imaginary frequencies.
        if sum(np.iscomplex(self.vib_energies)):
            raise ValueError('Imaginary vibrational energies are present.')
        else:
            self.vib_energies = np.real(self.vib_energies)  # clear +0.j 

        if electronicenergy:
            self.electronicenergy = electronicenergy
        else:
            self.electronicenergy = 0.

    def get_internal_energy(self, temperature, verbose=True):
        """Returns the internal energy, in eV, in the harmonic approximation
        at a specified temperature (K)."""

        self.verbose = verbose
        write = self._vprint
        fmt = '%-15s%13.3f eV'
        write('Internal energy components at T = %.2f K:' % temperature)
        write('='*31)

        U = 0.

        write(fmt % ('E_elec', self.electronicenergy))
        U += self.electronicenergy

        zpe = self.get_ZPE_correction()
        write(fmt % ('E_ZPE', zpe))
        U += zpe

        dU_v = self._vibrational_energy_contribution(temperature)
        write(fmt % ('Cv_harm (0->T)', dU_v))
        U += dU_v

        write('-'*31)
        write(fmt % ('U', U))
        write('='*31)
        return U

    def get_entropy(self, temperature, verbose=True):
        """Returns the entropy, in eV/K, in the harmonic approximation
        at a specified temperature (K)."""

        self.verbose = verbose
        write = self._vprint
        fmt = '%-15s%13.7f eV/K%13.3f eV'
        write('Entropy components at T = %.2f K:' % temperature)
        write('='*49)
        write('%15s%13s     %13s'%('', 'S', 'T*S'))

        S = 0.

        S_v = self._vibrational_entropy_contribution(temperature)
        write(fmt%('S_harm', S_v, S_v*temperature))
        S += S_v

        write('-'*49)
        write(fmt%('S', S, S*temperature))
        write('='*49)
        return S

    def get_free_energy(self, temperature, verbose=True):
        """Returns the free energy, in eV, in the harmonic approximation
        at a specified temperature (K)."""

        self.verbose = True
        write = self._vprint

        U = self.get_internal_energy(temperature, verbose=verbose)
        write('')
        S = self.get_entropy(temperature, verbose=verbose)
        G = U - temperature * S

        write('')
        write('Free energy components at T = %.2f K:' % temperature)
        write('='*23)
        fmt = '%5s%15.3f eV'
        write(fmt % ('U', U))
        write(fmt % ('-T*S', -temperature * S))
        write('-'*23)
        write(fmt % ('G', G))
        write('='*23)
        return G


class IdealGasThermo(ThermoChem):
    """Class for calculating thermodynamic properties of a molecule
    based on statistical mechanical treatments in the ideal gas
    approximation.

    Inputs for enthalpy calculations:

    vib_energies : list
        a list of the vibrational energies of the molecule (e.g., from
        ase.vibrations.Vibrations.get_energies). The number of vibrations
        used is automatically calculated by the geometry and the number of
        atoms. If more are specified than are needed, then the lowest 
        numbered vibrations are neglected. If either atoms or natoms is 
        unspecified, then uses the entire list. Units are eV.
    geometry : 'monatomic', 'linear', or 'nonlinear'
        geometry of the molecule
    electronicenergy : float
        the electronic energy in eV
        (If electronicenergy is unspecified, then the methods of this
        class can be interpreted as the enthalpy and free energy
        corrections.)
    natoms : integer
        the number of atoms, used along with 'geometry' to determine how
        many vibrations to use. (Not needed if an atoms object is supplied
        in 'atoms' or if the user desires the entire list of vibrations
        to be used.)

    Extra inputs needed for for entropy / free energy calculations:

    atoms : an ASE atoms object
        used to calculate rotational moments of inertia and molecular mass
    symmetrynumber : integer
        symmetry number of the molecule. See, for example, Table 10.1 and
        Appendix B of C. Cramer "Essentials of Computational Chemistry",
        2nd Ed.
    spin : float
        the total electronic spin. (0 for molecules in which all electrons
        are paired, 0.5 for a free radical with a single unpaired electron,
        1.0 for a triplet with two unpaired electrons, such as O_2.)
    """

    def __init__(self, vib_energies, geometry, electronicenergy=None,
                 atoms=None, symmetrynumber=None, spin=None, natoms=None):
        if electronicenergy == None:
            self.electronicenergy = 0.
        else:
            self.electronicenergy = electronicenergy
        self.geometry = geometry
        self.atoms = atoms
        self.sigma = symmetrynumber
        self.spin = spin
        if natoms == None:
            if atoms:
                natoms = len(atoms)
        # Cut the vibrations to those needed from the geometry.
        if natoms:
            if geometry == 'nonlinear':
                self.vib_energies = vib_energies[-(3*natoms-6):]
            elif geometry == 'linear':
                self.vib_energies = vib_energies[-(3*natoms-5):]
            elif geometry == 'monatomic':
                self.vib_energies = []
        else:
            self.vib_energies = vib_energies
        # Make sure no imaginary frequencies remain.
        if sum(np.iscomplex(self.vib_energies)):
            raise ValueError('Imaginary frequencies are present.')
        else:
            self.vib_energies = np.real(self.vib_energies) # clear +0.j 
        self.referencepressure = 101325. # Pa

    def get_enthalpy(self, temperature, verbose=True):
        """Returns the enthalpy, in eV, in the ideal gas approximation
        at a specified temperature (K)."""

        self.verbose = verbose
        write = self._vprint
        fmt = '%-15s%13.3f eV'
        write('Enthalpy components at T = %.2f K:' % temperature)
        write('='*31)

        H = 0.

        write(fmt % ('E_elec', self.electronicenergy))
        H += self.electronicenergy

        zpe = self.get_ZPE_correction()
        write(fmt % ('E_ZPE', zpe))
        H += zpe

        Cv_t = 3./2. * units.kB # translational heat capacity (3-d gas)
        write(fmt % ('Cv_trans (0->T)', Cv_t * temperature))
        H += Cv_t * temperature

        if self.geometry == 'nonlinear': # rotational heat capacity
            Cv_r = 3./2.*units.kB
        elif self.geometry == 'linear':
            Cv_r = units.kB
        elif self.geometry == 'monatomic':
            Cv_r = 0.
        write(fmt % ('Cv_rot (0->T)', Cv_r * temperature))
        H += Cv_r*temperature

        dH_v = self._vibrational_energy_contribution(temperature)
        write(fmt % ('Cv_vib (0->T)', dH_v))
        H += dH_v

        Cp_corr = units.kB * temperature
        write(fmt % ('(C_v -> C_p)', Cp_corr))
        H += Cp_corr

        write('-'*31)
        write(fmt % ('H', H))
        write('='*31)
        return H

    def get_entropy(self, temperature, pressure, verbose=True):
        """Returns the entropy, in eV/K, in the ideal gas approximation
        at a specified temperature (K) and pressure (Pa)."""

        if self.atoms == None or self.sigma == None or self.spin == None:
            raise RuntimeError('atoms, symmetrynumber, and spin must be '
                               'specified for entropy and free energy '
                               'calculations.')
        self.verbose = verbose
        write = self._vprint
        fmt = '%-15s%13.7f eV/K%13.3f eV'
        write('Entropy components at T = %.2f K and P = %.1f Pa:' %
               (temperature, pressure))
        write('=' * 49)
        write('%15s%13s     %13s' % ('', 'S', 'T*S'))

        S = 0.0

        # Translational entropy (term inside the log is in SI units).
        mass = sum(self.atoms.get_masses()) * units._amu  # kg/molecule
        S_t = (2 * np.pi * mass * units._k *
               temperature / units._hplanck**2)**(3.0 / 2)
        S_t *= units._k * temperature / self.referencepressure 
        S_t = units.kB * (np.log(S_t) + 5.0 / 2.0)
        write(fmt % ('S_trans (1 atm)', S_t, S_t * temperature))
        S += S_t

        # Rotational entropy (term inside the log is in SI units).
        if self.geometry == 'monatomic':
            S_r = 0.0
        elif self.geometry == 'nonlinear':
            inertias = (rotationalinertia(self.atoms) * units._amu /
                        (10.0**10)**2)  # kg m^2
            S_r = np.sqrt(np.pi * np.product(inertias)) / self.sigma
            S_r *= (8.0 * np.pi**2 * units._k * temperature /
                    units._hplanck**2)**(3.0 / 2.0)
            S_r = units.kB * (np.log(S_r) + 3.0 / 2.0)
        elif self.geometry == 'linear':
            inertias = (rotationalinertia(self.atoms) * units._amu /
                        (10.0**10)**2)  # kg m^2
            inertia = max(inertias)  # should be two identical and one zero
            S_r = (8 * np.pi**2 * inertia * units._k * temperature /
                   self.sigma / units._hplanck**2)
            S_r = units.kB * (np.log(S_r) + 1.)
        write(fmt % ('S_rot', S_r, S_r * temperature))
        S += S_r

        # Electronic entropy.
        S_e = units.kB * np.log(2 * self.spin + 1)
        write(fmt % ('S_elec', S_e, S_e * temperature))
        S += S_e

        # Vibrational entropy.
        S_v = self._vibrational_entropy_contribution(temperature)
        write(fmt % ('S_vib', S_v, S_v * temperature))
        S += S_v

        # Pressure correction to translational entropy.
        S_p = - units.kB * np.log(pressure / self.referencepressure)
        write(fmt % ('S (1 atm -> P)', S_p, S_p * temperature))
        S += S_p

        write('-' * 49)
        write(fmt % ('S', S, S * temperature))
        write('=' * 49)
        return S

    def get_free_energy(self, temperature, pressure, verbose=True):
        """Returns the free energy, in eV, in the ideal gas approximation
        at a specified temperature (K) and pressure (Pa)."""

        self.verbose = verbose
        write = self._vprint

        H = self.get_enthalpy(temperature, verbose=verbose)
        write('')
        S = self.get_entropy(temperature, pressure, verbose=verbose)
        G = H - temperature * S

        write('')
        write('Free energy components at T = %.2f K and P = %.1f Pa:' %
               (temperature, pressure))
        write('=' * 23)
        fmt = '%5s%15.3f eV'
        write(fmt % ('H', H))
        write(fmt % ('-T*S', -temperature * S))
        write('-' * 23)
        write(fmt % ('G', G))
        write('=' * 23)
        return G
