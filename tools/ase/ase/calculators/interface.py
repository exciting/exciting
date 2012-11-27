import numpy as np


class Calculator:
    """Class for demonstrating the ASE-calculator interface.

    A good implementation of a calculator should store a copy of the
    atoms object used for the last calculation.  When one of the
    *get_potential_energy*, *get_forces*, or *get_stress* methods is
    called, the calculator should check if anything has changed since
    the last calculation and only do the calculation if it's really
    needed.  The Atoms class implements the methods *__eq__* and
    *__ne__* that can be used for checking identity (using *==* and
    *!=*): Two sets of atoms are considered identical if they have the
    same positions, atomic numbers, unit cell and periodic boundary
    conditions."""

    def get_potential_energy(self, atoms=None, force_consistent=False):
        """Return total energy.
        
        Both the energy extrapolated to zero Kelvin and the energy
        consistent with the forces (the free energy) can be
        returned."""
        return 0.0
        
    def get_forces(self, atoms):
        """Return the forces."""
        return np.zeros((len(atoms), 3))
                        
    def get_stress(self, atoms):
        """Return the stress."""
        return np.zeros((3, 3))

    def calculation_required(self, atoms, quantities):
        """Check if a calculation is required.

        Check if the quantities in the *quantities* list have already
        been calculated for the atomic configuration *atoms*.  The
        quantities can be one or more of: 'energy', 'forces', 'stress',
        and 'magmoms'."""
        return False

    def set_atoms(self, atoms):
        """Let the calculator know the atoms.

        This method is optional.  If it exists, it will be called when
        the *Atoms.set_calculator()* method is called.

        *Don't* store a reference to *atoms* - that will create a
         cyclic reference loop!  Store a copy instead."""
        self.atoms = atoms.copy()


class DFTCalculator(Calculator):
    """Class for demonstrating the ASE interface to DFT-calculators."""
    
    def get_number_of_bands(self):
        """Return the number of bands."""
        return 42
  
    def get_xc_functional(self):
        """Return the XC-functional identifier.
        
        'LDA', 'PBE', ..."""
        return 'LDA'
 
    def get_bz_k_points(self):
        """Return all the k-points in the 1. Brillouin zone.

        The coordinates are relative to reciprocal latice vectors."""
        return np.zeros((1, 3))
 
    def get_number_of_spins(self):
        """Return the number of spins in the calculation.

        Spin-paired calculations: 1, spin-polarized calculation: 2."""
        return 1

    def get_spin_polarized(self):
        """Is it a spin-polarized calculation?"""
        return False
    
    def get_ibz_k_points(self):
        """Return k-points in the irreducible part of the Brillouin zone.

        The coordinates are relative to reciprocal latice vectors."""
        return np.zeros((1, 3))

    def get_k_point_weights(self):
        """Weights of the k-points. 
        
        The sum of all weights is one."""
        return np.ones(1)

    def get_pseudo_density(self, spin=None, pad=True):
        """Return pseudo-density array.
        
        If *spin* is not given, then the total density is returned.
        Otherwise, the spin up or down density is returned (spin=0 or
        1)."""
        return np.zeros((40, 40, 40))

    def get_effective_potential(self, spin=0, pad=True):
        """Return pseudo-effective-potential array."""
        return np.zeros((40, 40, 40))

    def get_pseudo_wave_function(self, band=0, kpt=0, spin=0, broadcast=True,
                                 pad=True):
        """Return pseudo-wave-function array."""
        return np.zeros((40, 40, 40))

    def get_eigenvalues(self, kpt=0, spin=0):
        """Return eigenvalue array."""
        return np.arange(42, float)
    
    def get_occupation_numbers(self, kpt=0, spin=0):
        """Return occupation number array."""
        return np.ones(42)
        
    def get_fermi_level(self):
        """Return the Fermi level."""
        return 0.0

    def initial_wannier(self, initialwannier, kpointgrid, fixedstates,
                        edf, spin, nbands):
        """Initial guess for the shape of wannier functions.

        Use initial guess for wannier orbitals to determine rotation
        matrices U and C.
        """
        raise NotImplementedError

        return c, U

    def get_wannier_localization_matrix(self, nbands, dirG, kpoint,
                                        nextkpoint, G_I, spin):
        """Calculate integrals for maximally localized Wannier functions."""

    def get_magnetic_moment(self, atoms=None):
        """Return the total magnetic moment."""
        return self.occupation.magmom

    def get_number_of_grid_points(self):
        """Return the shape of arrays."""
        return self.gd.N_c


