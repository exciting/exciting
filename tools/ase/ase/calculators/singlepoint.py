import numpy as np


class SinglePointCalculator:
    """Special calculator for a single configuration.

    Used to remember the energy, force and stress for a given
    configuration.  If the positions, atomic numbers, unit cell, or
    boundary conditions are changed, then asking for
    energy/forces/stress will raise an exception."""
    
    def __init__(self, energy, forces, stress, magmoms, atoms):
        """Save energy, forces and stresses for the current configuration."""
        self.energy = energy
        if forces is not None:
            forces = np.array(forces, float)
        self.forces = forces
        if stress is not None:
            stress = np.array(stress, float)
        self.stress = stress
        if magmoms is not None:
            magmoms = np.array(magmoms, float)
        self.magmoms = magmoms
        self.atoms = atoms.copy()

    def calculation_required(self, atoms, quantities):
        ok = self.atoms == atoms
        return ('forces' in quantities and (self.forces is None or not ok) or
                'energy' in quantities and (self.energy is None or not ok) or
                'stress' in quantities and (self.stress is None or not ok) or
                'magmoms' in quantities and (self.magmoms is None or not ok))

    def update(self, atoms):
        if self.atoms != atoms:
            raise RuntimeError('Energy, forces and stress no longer correct.')

    def get_potential_energy(self, atoms=None):
        if atoms is not None:
            self.update(atoms)
        if self.energy is None:
            raise RuntimeError('No energy.')
        return self.energy

    def get_forces(self, atoms):
        self.update(atoms)
        if self.forces is None:
            raise RuntimeError('No forces.')
        return self.forces

    def get_stress(self, atoms):
        self.update(atoms)
        if self.stress is None:
            raise NotImplementedError
        return self.stress

    def get_spin_polarized(self):
        return self.magmoms is not None and self.magmoms.any()

    def get_magnetic_moments(self, atoms=None):
        if atoms is not None:
            self.update(atoms)
        if self.magmoms is not None:
            return self.magmoms
        else:
            return np.zeros(len(self.positions))

    def get_dipole_moment(self):
        if hasattr(self, 'dipole'):
            return self.dipole
        return None

    def set_dipole_moment(self, dipole):
        self.dipole = dipole

class SinglePointKPoint:
    def __init__(self, kpt, spin):
        self.k = kpt
        self.s = spin
        self.eps_n = []
        self.f_n = []

class SinglePointDFTCalculator(SinglePointCalculator):
    def __init__(self, energy, forces, stress, magmoms, atoms,
                 eFermi=None):
        SinglePointCalculator.__init__(self, energy, forces, stress, 
                                       magmoms, atoms)
        if eFermi is not None:
            self.eFermi = eFermi
        self.kpts = None

    def get_fermi_level(self):
        """Return the Fermi-level(s)."""
        return self.eFermi

    def get_bz_k_points(self):
        """Return the k-points."""
        if self.kpts is not None:
            # we assume that only the gamma point is defined
            return np.zeros((1, 3))
        return None

    def get_number_of_spins(self):
        """Return the number of spins in the calculation.

        Spin-paired calculations: 1, spin-polarized calculation: 2."""
        if self.kpts is not None:
            # we assume that only the gamma point is defined
            return len(self.kpts)
        return None

    def get_spin_polarized(self):
        """Is it a spin-polarized calculation?"""
        nos = self.get_number_of_spins()
        if nos is not None:
            return nos == 2
        return None
    
    def get_ibz_k_points(self):
        """Return k-points in the irreducible part of the Brillouin zone."""
        return self.get_bz_k_points()

    def get_occupation_numbers(self, kpt=0, spin=0):
        """Return occupation number array."""
        # we assume that only the gamma point is defined
        assert(kpt == 0)
        if self.kpts is not None:
            for kpt in self.kpts:
                if kpt.s == spin:
                    return kpt.f_n
        return None

    def get_eigenvalues(self, kpt=0, spin=0):
        """Return eigenvalue array."""
        # we assume that only the gamma point is defined
        assert(kpt == 0)
        if self.kpts is not None:
            for kpt in self.kpts:
                if kpt.s == spin:
                    return kpt.eps_n
        return None
