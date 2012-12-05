"""Module for calculating phonons of periodic systems."""

import sys
import pickle
from math import sin, pi, sqrt
from os import remove
from os.path import isfile

import numpy as np
import numpy.linalg as la
import numpy.fft as fft

has_spglib = False
try:
    from pyspglib import spglib
    has_spglib = True
except ImportError:
    pass

import ase.units as units
from ase.parallel import rank, barrier
from ase.dft import monkhorst_pack
from ase.io.trajectory import PickleTrajectory

class Displacement:
    """Abstract base class for phonon and el-ph supercell calculations.

    Both phonons and the electron-phonon interaction in periodic systems can be
    calculated with the so-called finite-displacement method where the
    derivatives of the total energy and effective potential are obtained from
    finite-difference approximations, i.e. by displacing the atoms. This class
    provides the required functionality for carrying out the calculations for
    the different displacements in its ``run`` member function.

    Derived classes must overwrite the ``__call__`` member function which is
    called for each atomic displacement.
    
    """

    def __init__(self, atoms, calc=None, supercell=(1, 1, 1), name=None,
                 delta=0.01, refcell=None):
        """Init with an instance of class ``Atoms`` and a calculator.

        Parameters
        ----------
        atoms: Atoms object
            The atoms to work on.
        calc: Calculator
            Calculator for the supercell calculation.
        supercell: tuple
            Size of supercell given by the number of repetitions (l, m, n) of
            the small unit cell in each direction.
        name: str
            Base name to use for files.
        delta: float
            Magnitude of displacement in Ang.
        refcell: str
            Reference cell in which the atoms will be displaced. If ``None``,
            corner cell in supercell is used. If ``str``, cell in the center of
            the supercell is used.

        """

        # Store atoms and calculator
        self.atoms = atoms
        self.calc = calc
        
        # Displace all atoms in the unit cell by default
        self.indices = range(len(atoms))
        self.name = name
        self.delta = delta
        self.N_c = supercell

        # Reference cell offset
        if refcell is None:
            # Corner cell
            self.offset = 0
        else:
            # Center cell
            N_c = self.N_c
            self.offset = N_c[0] // 2 * (N_c[1] * N_c[2]) + N_c[1] // \
                          2 * N_c[2] + N_c[2] // 2        

    def __call__(self, *args, **kwargs):
        """Member function called in the ``run`` function."""

        raise NotImplementedError("Implement in derived classes!.")
    
    def set_atoms(self, atoms):
        """Set the atoms to vibrate.

        Parameters
        ----------
        atoms: list
            Can be either a list of strings, ints or ...
            
        """
        
        assert isinstance(atoms, list)
        assert len(atoms) <= len(self.atoms)
        
        if isinstance(atoms[0], str):
            assert np.all([isinstance(atom, str) for atom in atoms])
            sym_a = self.atoms.get_chemical_symbols()
            # List for atomic indices
            indices = []
            for type in atoms:
                indices.extend([a for a, atom in enumerate(sym_a)
                                if atom == type])
        else:
            assert np.all([isinstance(atom, int) for atom in atoms])
            indices = atoms

        self.indices = indices

    def lattice_vectors(self):
        """Return lattice vectors for cells in the supercell."""

        # Lattice vectors relevative to the reference cell
        R_cN = np.indices(self.N_c).reshape(3, -1)
        N_c = np.array(self.N_c)[:, np.newaxis]
        if self.offset == 0:
            R_cN += N_c // 2
            R_cN %= N_c
        R_cN -= N_c // 2

        return R_cN
    
    def run(self):
        """Run the calculations for the required displacements.

        This will do a calculation for 6 displacements per atom, +-x, +-y, and
        +-z. Only those calculations that are not already done will be
        started. Be aware that an interrupted calculation may produce an empty
        file (ending with .pckl), which must be deleted before restarting the
        job. Otherwise the calculation for that displacement will not be done.

        """

        # Atoms in the supercell -- repeated in the lattice vector directions
        # beginning with the last
        atoms_N = self.atoms * self.N_c
        
        # Set calculator if provided
        assert self.calc is not None, "Provide calculator in __init__ method"
        atoms_N.set_calculator(self.calc)
        
        # Do calculation on equilibrium structure
        filename = self.name + '.eq.pckl'
        
        if not isfile(filename):
            # Wait for all ranks to enter
            barrier()
            # Create file
            if rank == 0:
                fd = open(filename, 'w')
                fd.close()

            # Call derived class implementation of __call__
            output = self.__call__(atoms_N)
            # Write output to file
            if rank == 0:
                fd = open(filename, 'w')
                pickle.dump(output, fd)
                sys.stdout.write('Writing %s\n' % filename)
                fd.close()
            sys.stdout.flush()

        # Positions of atoms to be displaced in the reference cell
        natoms = len(self.atoms)
        offset = natoms * self.offset
        pos = atoms_N.positions[offset: offset + natoms].copy()
        
        # Loop over all displacements
        for a in self.indices:
            for i in range(3):
                for sign in [-1, 1]:
                    # Filename for atomic displacement
                    filename = '%s.%d%s%s.pckl' % \
                               (self.name, a, 'xyz'[i], ' +-'[sign])
                    # Wait for ranks before checking for file
                    # barrier()                    
                    if isfile(filename):
                        # Skip if already done
                        continue
                    # Wait for ranks
                    barrier()
                    if rank == 0:
                        fd = open(filename, 'w')
                        fd.close()

                    # Update atomic positions
                    atoms_N.positions[offset + a, i] = \
                        pos[a, i] + sign * self.delta
                    
                    # Call derived class implementation of __call__
                    output = self.__call__(atoms_N)
                    # Write output to file    
                    if rank == 0:
                        fd = open(filename, 'w')
                        pickle.dump(output, fd)
                        sys.stdout.write('Writing %s\n' % filename)
                        fd.close()
                    sys.stdout.flush()
                    # Return to initial positions
                    atoms_N.positions[offset + a, i] = pos[a, i]
        
    def clean(self):
        """Delete generated pickle files."""
        
        if isfile(self.name + '.eq.pckl'):
            remove(self.name + '.eq.pckl')
        
        for a in self.indices:
            for i in 'xyz':
                for sign in '-+':
                    name = '%s.%d%s%s.pckl' % (self.name, a, i, sign)
                    if isfile(name):
                        remove(name)


class Phonons(Displacement):
    """Class for calculating phonon modes using the finite displacement method.

    The matrix of force constants is calculated from the finite difference
    approximation to the first-order derivative of the atomic forces as::
    
                            2             nbj   nbj
                nbj        d E           F-  - F+
               C     = ------------ ~  -------------  ,
                mai     dR   dR          2 * delta
                          mai  nbj       

    where F+/F- denotes the force in direction j on atom nb when atom ma is
    displaced in direction +i/-i. The force constants are related by various
    symmetry relations. From the definition of the force constants it must
    be symmetric in the three indices mai::

                nbj    mai         bj        ai
               C    = C      ->   C  (R ) = C  (-R )  .
                mai    nbj         ai  n     bj   n

    As the force constants can only depend on the difference between the m and
    n indices, this symmetry is more conveniently expressed as shown on the
    right hand-side.

    The acoustic sum-rule::

                           _ _
                aj         \    bj    
               C  (R ) = -  )  C  (R )
                ai  0      /__  ai  m
                          (m, b)
                            !=
                          (0, a)
                        
    Ordering of the unit cells illustrated here for a 1-dimensional system (in
    case ``refcell=None`` in constructor!):
    
    ::
    
               m = 0        m = 1        m = -2        m = -1
           -----------------------------------------------------
           |            |            |            |            |
           |        * b |        *   |        *   |        *   |
           |            |            |            |            |
           |   * a      |   *        |   *        |   *        |
           |            |            |            |            |
           -----------------------------------------------------
       
    Example:

    >>> from ase.structure import bulk
    >>> from ase.phonons import Phonons
    >>> from gpaw import GPAW, FermiDirac
    >>> atoms = bulk('Si', 'diamond', a=5.4)
    >>> calc = GPAW(kpts=(5, 5, 5),
                    h=0.2,
                    occupations=FermiDirac(0.))
    >>> ph = Phonons(atoms, calc, supercell=(5, 5, 5))
    >>> ph.run()
    >>> ph.read(method='frederiksen', acoustic=True)

    """

    def __init__(self, *args, **kwargs):
        """Initialize with base class args and kwargs."""

        if 'name' not in kwargs.keys():
            kwargs['name'] = "phonon"
            
        Displacement.__init__(self, *args, **kwargs)
        
        # Attributes for force constants and dynamical matrix in real space
        self.C_N = None  # in units of eV / Ang**2 
        self.D_N = None  # in units of eV / Ang**2 / amu
        
        # Attributes for born charges and static dielectric tensor
        self.Z_avv = None
        self.eps_vv = None
   
    def __call__(self, atoms_N):
        """Calculate forces on atoms in supercell."""

        # Calculate forces
        forces = atoms_N.get_forces()

        return forces

    def check_eq_forces(self):
        """Check maximum size of forces in the equilibrium structure."""

        fname = '%s.eq.pckl' % self.name
        feq_av = pickle.load(open(fname))

        fmin = feq_av.max()
        fmax = feq_av.min()
        i_min = np.where(feq_av == fmin)
        i_max = np.where(feq_av == fmax)

        return fmin, fmax, i_min, i_max
    
    def read_born_charges(self, name=None, neutrality=True):
        """Read Born charges and dieletric tensor from pickle file.

        The charge neutrality sum-rule::
    
                   _ _
                   \    a    
                    )  Z   = 0
                   /__  ij
                    a
                              
        Parameters
        ----------
        neutrality: bool
            Restore charge neutrality condition on calculated Born effective
            charges. 

        """

        # Load file with Born charges and dielectric tensor for atoms in the
        # unit cell
        if name is None:
            filename = '%s.born.pckl' % self.name
        else:
            filename = name
            
        fd = open(filename)
        Z_avv, eps_vv = pickle.load(fd)
        fd.close()

        # Neutrality sum-rule
        if neutrality:
            Z_mean = Z_avv.sum(0) / len(Z_avv)
            Z_avv -= Z_mean
        
        self.Z_avv = Z_avv[self.indices]
        self.eps_vv = eps_vv
        
    def read(self, method='Frederiksen', symmetrize=3, acoustic=True,
             cutoff=None, born=False, **kwargs):
        """Read forces from pickle files and calculate force constants.

        Extra keyword arguments will be passed to ``read_born_charges``.
        
        Parameters
        ----------
        method: str
            Specify method for evaluating the atomic forces.
        symmetrize: int
            Symmetrize force constants (see doc string at top) when
            ``symmetrize != 0`` (default: 3). Since restoring the acoustic sum
            rule breaks the symmetry, the symmetrization must be repeated a few
            times until the changes a insignificant. The integer gives the
            number of iterations that will be carried out.
        acoustic: bool
            Restore the acoustic sum rule on the force constants.
        cutoff: None or float
            Zero elements in the dynamical matrix between atoms with an
            interatomic distance larger than the cutoff.
        born: bool
            Read in Born effective charge tensor and high-frequency static
            dielelctric tensor from file.
            
        """

        method = method.lower()
        assert method in ['standard', 'frederiksen']
        if cutoff is not None:
            cutoff = float(cutoff)
            
        # Read Born effective charges and optical dielectric tensor
        if born:
            self.read_born_charges(**kwargs)
        
        # Number of atoms
        natoms = len(self.indices)
        # Number of unit cells
        N = np.prod(self.N_c)
        # Matrix of force constants as a function of unit cell index in units
        # of eV / Ang**2
        C_xNav = np.empty((natoms * 3, N, natoms, 3), dtype=float)

        # Loop over all atomic displacements and calculate force constants
        for i, a in enumerate(self.indices):
            for j, v in enumerate('xyz'):
                # Atomic forces for a displacement of atom a in direction v
                basename = '%s.%d%s' % (self.name, a, v)
                fminus_av = pickle.load(open(basename + '-.pckl'))
                fplus_av = pickle.load(open(basename + '+.pckl'))
                
                if method == 'frederiksen':
                    fminus_av[a] -= fminus_av.sum(0)
                    fplus_av[a]  -= fplus_av.sum(0)
                    
                # Finite difference derivative
                C_av = fminus_av - fplus_av
                C_av /= 2 * self.delta

                # Slice out included atoms
                C_Nav = C_av.reshape((N, len(self.atoms), 3))[:, self.indices]
                index = 3*i + j                
                C_xNav[index] = C_Nav

        # Make unitcell index the first and reshape
        C_N = C_xNav.swapaxes(0 ,1). reshape((N,) + (3 * natoms, 3 * natoms))

        # Cut off before symmetry and acoustic sum rule are imposed
        if cutoff is not None:
            self.apply_cutoff(C_N, cutoff)
            
        # Symmetrize force constants
        if symmetrize:
            for i in range(symmetrize):
                # Symmetrize
                C_N = self.symmetrize(C_N)
                # Restore acoustic sum-rule
                if acoustic:
                    self.acoustic(C_N)
                else:
                    break
             
        # Store force constants and dynamical matrix
        self.C_N = C_N
        self.D_N = C_N.copy()
        
        # Add mass prefactor
        m_a = self.atoms.get_masses()
        self.m_inv_x = np.repeat(m_a[self.indices]**-0.5, 3)
        M_inv = np.outer(self.m_inv_x, self.m_inv_x)
        for D in self.D_N:
            D *= M_inv

    def symmetrize(self, C_N):
        """Symmetrize force constant matrix."""

        # Number of atoms
        natoms = len(self.indices)
        # Number of unit cells
        N = np.prod(self.N_c)

        # Reshape force constants to (l, m, n) cell indices
        C_lmn = C_N.reshape(self.N_c + (3 * natoms, 3 * natoms))

        # Shift reference cell to center index 
        if self.offset == 0:
            C_lmn = fft.fftshift(C_lmn, axes=(0, 1, 2)).copy()
        # Make force constants symmetric in indices -- in case of an even
        # number of unit cells don't include the first
        i, j, k = np.asarray(self.N_c) % 2 - 1
        C_lmn[i:, j:, k:] *= 0.5
        C_lmn[i:, j:, k:] += \
                  C_lmn[i:, j:, k:][::-1, ::-1, ::-1].transpose(0, 1, 2, 4, 3).copy()
        if self.offset == 0:
            C_lmn = fft.ifftshift(C_lmn, axes=(0, 1, 2)).copy()

        # Change to single unit cell index shape
        C_N = C_lmn.reshape((N, 3 * natoms, 3 * natoms))

        return C_N
       
    def acoustic(self, C_N):
        """Restore acoustic sumrule on force constants."""

        # Number of atoms
        natoms = len(self.indices)
        # Copy force constants
        C_N_temp = C_N.copy()
        
        # Correct atomic diagonals of R_m = (0, 0, 0) matrix
        for C in C_N_temp:
            for a in range(natoms):
                for a_ in range(natoms):
                    C_N[self.offset, 3*a: 3*a + 3, 3*a: 3*a + 3] -= \
                                     C[3*a: 3*a+3, 3*a_: 3*a_+3]
                    
    def apply_cutoff(self, D_N, r_c):
        """Zero elements for interatomic distances larger than the cutoff.

        Parameters
        ----------
        D_N: ndarray
            Dynamical/force constant matrix.
        r_c: float
            Cutoff in Angstrom.

        """

        # Number of atoms and primitive cells
        natoms = len(self.indices)
        N = np.prod(self.N_c)
        # Lattice vectors
        R_cN = self.lattice_vectors()
        # Reshape matrix to individual atomic and cartesian dimensions
        D_Navav = D_N.reshape((N, natoms, 3, natoms, 3))

        # Cell vectors
        cell_vc = self.atoms.cell.transpose()
        # Atomic positions in reference cell
        pos_av = self.atoms.get_positions()
        
        # Zero elements with a distance to atoms in the reference cell
        # larger than the cutoff
        for n in range(N):
            # Lattice vector to cell
            R_v = np.dot(cell_vc, R_cN[:, n])
            # Atomic positions in cell
            posn_av = pos_av + R_v
            # Loop over atoms and zero elements
            for i, a in enumerate(self.indices):
                dist_a = np.sqrt(np.sum((pos_av[a] - posn_av)**2, axis=-1))
                # Atoms where the distance is larger than the cufoff
                i_a = dist_a > r_c #np.where(dist_a > r_c)
                # Zero elements
                D_Navav[n, i, :, i_a, :] = 0.0
            # print ""
            
    def get_force_constant(self):
        """Return matrix of force constants."""

        assert self.C_N is not None
        
        return self.C_N
    
    def band_structure(self, path_kc, modes=False, born=False, verbose=True):
        """Calculate phonon dispersion along a path in the Brillouin zone.

        The dynamical matrix at arbitrary q-vectors is obtained by Fourier
        transforming the real-space force constants. In case of negative
        eigenvalues (squared frequency), the corresponding negative frequency
        is returned.

        Eigenvalues and modes are in units of eV and Ang/sqrt(amu),
        respectively.

        Parameters
        ----------
        path_kc: ndarray
            List of k-point coordinates (in units of the reciprocal lattice
            vectors) specifying the path in the Brillouin zone for which the
            dynamical matrix will be calculated.
        modes: bool
            Returns both frequencies and modes when True.
        born: bool
            Include non-analytic part given by the Born effective charges and
            the static part of the high-frequency dielectric tensor. This
            contribution to the force constant accounts for the splitting
            between the LO and TO branches for q -> 0.
        verbose: bool
            Print warnings when imaginary frequncies are detected.
        
        """

        assert self.D_N is not None
        if born:
            assert self.Z_avv is not None
            assert self.eps_vv is not None

        # Lattice vectors -- ordered as illustrated in class docstring
        R_cN = self.lattice_vectors()

        # Dynamical matrix in real-space
        D_N = self.D_N
        
        # Lists for frequencies and modes along path
        omega_kl = []
        u_kl = []

        # Reciprocal basis vectors for use in non-analytic contribution
        reci_vc = 2 * pi * la.inv(self.atoms.cell)
        # Unit cell volume in Bohr^3
        vol = abs(la.det(self.atoms.cell)) / units.Bohr**3

        for q_c in path_kc:
           
            # Add non-analytic part
            if born:
                # q-vector in cartesian coordinates
                q_v = np.dot(reci_vc, q_c)
                # Non-analytic contribution to force constants in atomic units
                qdotZ_av = np.dot(q_v, self.Z_avv).ravel()
                C_na = 4 * pi * np.outer(qdotZ_av, qdotZ_av) / \
                       np.dot(q_v, np.dot(self.eps_vv, q_v)) / vol
                self.C_na = C_na / units.Bohr**2 * units.Hartree                
                # Add mass prefactor and convert to eV / (Ang^2 * amu)
                M_inv = np.outer(self.m_inv_x, self.m_inv_x)                
                D_na = C_na * M_inv / units.Bohr**2 * units.Hartree
                self.D_na = D_na
                D_N = self.D_N + D_na / np.prod(self.N_c) 

            ## if np.prod(self.N_c) == 1:
            ## 
            ##     q_av = np.tile(q_v, len(self.indices))
            ##     q_xx = np.vstack([q_av]*len(self.indices)*3)
            ##     D_m += q_xx

            # Evaluate fourier sum
            phase_N = np.exp(-2.j * pi * np.dot(q_c, R_cN))
            D_q = np.sum(phase_N[:, np.newaxis, np.newaxis] * D_N, axis=0)

            if modes:
                omega2_l, u_xl = la.eigh(D_q, UPLO='U')
                # Sort eigenmodes according to eigenvalues (see below) and 
                # multiply with mass prefactor
                u_lx = (self.m_inv_x[:, np.newaxis] * 
                        u_xl[:, omega2_l.argsort()]).T.copy()
                u_kl.append(u_lx.reshape((-1, len(self.indices), 3)))
            else:
                omega2_l = la.eigvalsh(D_q, UPLO='U')

            # Sort eigenvalues in increasing order
            omega2_l.sort()
            # Use dtype=complex to handle negative eigenvalues
            omega_l = np.sqrt(omega2_l.astype(complex))

            # Take care of imaginary frequencies
            if not np.all(omega2_l >= 0.):
                indices = np.where(omega2_l < 0)[0]

                if verbose:
                    print ("WARNING, %i imaginary frequencies at "
                           "q = (% 5.2f, % 5.2f, % 5.2f) ; (omega_q =% 5.3e*i)"
                           % (len(indices), q_c[0], q_c[1], q_c[2],
                              omega_l[indices][0].imag))
                
                omega_l[indices] = -1 * np.sqrt(np.abs(omega2_l[indices].real))

            omega_kl.append(omega_l.real)

        # Conversion factor: sqrt(eV / Ang^2 / amu) -> eV
        s = units._hbar * 1e10 / sqrt(units._e * units._amu)
        omega_kl = s * np.asarray(omega_kl)
        
        if modes:
            return omega_kl, np.asarray(u_kl)
        
        return omega_kl

    def dos(self, kpts=(10, 10, 10), npts=1000, delta=1e-3, indices=None):
        """Calculate phonon dos as a function of energy.

        Parameters
        ----------
        qpts: tuple
            Shape of Monkhorst-Pack grid for sampling the Brillouin zone.
        npts: int
            Number of energy points.
        delta: float
            Broadening of Lorentzian line-shape in eV.
        indices: list
            If indices is not None, the atomic-partial dos for the specified
            atoms will be calculated.
            
        """

        # Monkhorst-Pack grid
        kpts_kc = monkhorst_pack(kpts)
        N = np.prod(kpts)
        # Get frequencies
        omega_kl = self.band_structure(kpts_kc)
        # Energy axis and dos
        omega_e = np.linspace(0., np.amax(omega_kl) + 5e-3, num=npts)
        dos_e = np.zeros_like(omega_e)
       
        # Sum up contribution from all q-points and branches
        for omega_l in omega_kl:
            diff_el = (omega_e[:, np.newaxis] - omega_l[np.newaxis, :])**2
            dos_el = 1. / (diff_el + (0.5 * delta)**2)
            dos_e += dos_el.sum(axis=1)

        dos_e *= 1./(N * pi) * 0.5*delta
        
        return omega_e, dos_e
    
    def write_modes(self, q_c, branches=0, kT=units.kB*300, born=False,
                    repeat=(1, 1, 1), nimages=30, center=False):
        """Write modes to trajectory file.

        Parameters
        ----------
        q_c: ndarray
            q-vector of modes.
        branches: int or list
            Branch index of modes.
        kT: float
            Temperature in units of eV. Determines the amplitude of the atomic
            displacements in the modes.
        born: bool
            Include non-analytic contribution to the force constants at q -> 0.
        repeat: tuple
            Repeat atoms (l, m, n) times in the directions of the lattice
            vectors. Displacements of atoms in repeated cells carry a Bloch
            phase factor given by the q-vector and the cell lattice vector R_m.
        nimages: int
            Number of images in an oscillation.
        center: bool
            Center atoms in unit cell if True (default: False).
            
        """

        if isinstance(branches, int):
            branch_l = [branches]
        else:
            branch_l = list(branches)

        # Calculate modes
        omega_l, u_l = self.band_structure([q_c], modes=True, born=born)
        # Repeat atoms
        atoms = self.atoms * repeat
        # Center
        if center:
            atoms.center()
        
        # Here ``Na`` refers to a composite unit cell/atom dimension
        pos_Nav = atoms.get_positions()
        # Total number of unit cells
        N = np.prod(repeat)

        # Corresponding lattice vectors R_m
        R_cN = np.indices(repeat).reshape(3, -1)
        # Bloch phase
        phase_N = np.exp(2.j * pi * np.dot(q_c, R_cN))
        phase_Na = phase_N.repeat(len(self.atoms))

        for l in branch_l:
            
            omega = omega_l[0, l]
            u_av = u_l[0, l]

            # Mean displacement of a classical oscillator at temperature T
            u_av *= sqrt(kT) / abs(omega)

            mode_av = np.zeros((len(self.atoms), 3), dtype=complex)
            # Insert slice with atomic displacements for the included atoms
            mode_av[self.indices] = u_av
            # Repeat and multiply by Bloch phase factor
            mode_Nav = (np.vstack(N * [mode_av]) * phase_Na[:, np.newaxis]).real
            
            traj = PickleTrajectory('%s.mode.%d.traj' % (self.name, l), 'w')
            
            for x in np.linspace(0, 2*pi, nimages, endpoint=False):
                atoms.set_positions(pos_Nav + sin(x) * mode_Nav)
                traj.write(atoms)
                
            traj.close()
