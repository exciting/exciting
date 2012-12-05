'''Constant pressure/stress and temperature dynamics.

Combined Nose-Hoover and Parrinello-Rahman dynamics, creating an NPT
(or N,stress,T) ensemble.

The method is the one proposed by Melchionna et al. [1] and later
modified by Melchionna [2].  The differential equations are integrated
using a centered difference method [3].

 1. S. Melchionna, G. Ciccotti and B. L. Holian, "Hoover NPT dynamics
    for systems varying in shape and size", Molecular Physics 78, p. 533
    (1993).

 2. S. Melchionna, "Constrained systems and statistical distribution",
    Physical Review E, 61, p. 6165 (2000).

 3. B. L. Holian, A. J. De Groot, W. G. Hoover, and C. G. Hoover,
    "Time-reversible equilibrium and nonequilibrium isothermal-isobaric
    simulations with centered-difference Stoermer algorithms.", Physical
    Review A, 41, p. 4552 (1990).
'''

__docformat__ = 'reStructuredText'

from numpy import *
import sys
import weakref
from ase.md.md import MolecularDynamics
#from ASE.Trajectories.NetCDFTrajectory import NetCDFTrajectory

# Delayed imports:  If the trajectory object is reading a special ASAP version
# of HooverNPT, that class is imported from Asap.Dynamics.NPTDynamics.

class NPT(MolecularDynamics):
    '''Constant pressure/stress and temperature dynamics.

    Combined Nose-Hoover and Parrinello-Rahman dynamics, creating an
    NPT (or N,stress,T) ensemble.

    The method is the one proposed by Melchionna et al. [1] and later
    modified by Melchionna [2].  The differential equations are integrated
    using a centered difference method [3].  See also NPTdynamics.tex

    The dynamics object is called with the following parameters:

    atoms
        The list of atoms.

    dt
        The timestep in units matching eV, A, u.

    temperature
        The desired temperature in eV.

    externalstress
        The external stress in eV/A^3.  Either a symmetric
        3x3 tensor, a 6-vector representing the same, or a
        scalar representing the pressure.  Note that the
        stress is positive in tension whereas the pressure is
        positive in compression: giving a scalar p is
        equivalent to giving the tensor (-p, -p, -p, 0, 0, 0).

    ttime
        Characteristic timescale of the thermostat.
        Set to None to disable the thermostat.

    pfactor
        A constant in the barostat differential equation.  If
        a characteristic barostat timescale of ptime is
        desired, set pfactor to ptime^2 * B (where B is the
        Bulk Modulus).  Set to None to disable the barostat.
        Typical metallic bulk moduli are of the order of
        100 GPa or 0.6 eV/A^3.

    mask=None
        Optional argument.  A tuple of three integers (0 or 1),
        indicating if the system can change size along the
        three Cartesian axes.  Set to (1,1,1) or None to allow
        a fully flexible computational box.  Set to (1,1,0)
        to disallow elongations along the z-axis etc.

    Useful parameter values:

    * The same timestep can be used as in Verlet dynamics, i.e. 5 fs is fine
      for bulk copper.

    * The ttime and pfactor are quite critical[4], too small values may
      cause instabilites and/or wrong fluctuations in T / p.  Too
      large values cause an oscillation which is slow to die.  Good
      values for the characteristic times seem to be 25 fs for ttime,
      and 75 fs for ptime (used to calculate pfactor), at least for
      bulk copper with 15000-200000 atoms.  But this is not well
      tested, it is IMPORTANT to monitor the temperature and
      stress/pressure fluctuations.
    
    It has the following methods:

    __call__(n)
        Perform n timesteps.
    initialize()
        Estimates the dynamic variables for time=-1 to start
        the algorithm.   This is automatically called before
        the first timestep.
    set_stress()
        Set the external stress.  Use with care.  It is
        preferable to set the right value when creating the
        object.
    set_mask()
        Change the mask.  Use with care, as you may "freeze"
        a fluctuation in the strain rate.
    get_gibbs_free_energy()
        Gibbs free energy is supposed to be preserved by this
        dynamics.  This is mainly intended as a diagnostic
        tool.
    
    References:
    
    1) S. Melchionna, G. Ciccotti and B. L. Holian, Molecular
       Physics 78, p. 533 (1993).

    2) S. Melchionna, Physical
       Review E 61, p. 6165 (2000).

    3) B. L. Holian, A. J. De Groot, W. G. Hoover, and C. G. Hoover,
       Physical Review A 41, p. 4552 (1990).

    4) F. D. Di Tolla and M. Ronchetti, Physical
       Review E 48, p. 1726 (1993).
    
    '''

    classname = "NPT"  # Used by the trajectory.
    def __init__(self, atoms, 
                 timestep, temperature, externalstress, ttime, pfactor,
                 mask=None, trajectory=None, logfile=None, loginterval=1):
        MolecularDynamics.__init__(self, atoms, timestep, trajectory,
                                   logfile, loginterval)
        #self.atoms = atoms
        #self.timestep = timestep
        self.zero_center_of_mass_momentum(verbose=1)
        self.temperature = temperature
        self.set_stress(externalstress)
        self.set_mask(mask)
        self.eta = zeros((3,3), float)
        self.zeta = 0.0
        self.zeta_integrated = 0.0
        self.initialized = 0
        self.ttime = ttime
        self.pfactor_given = pfactor
        self._calculateconstants()
        self.timeelapsed = 0.0
        self.frac_traceless = 1

    def set_temperature(self, temperature):
        self.temperature = temperature
        self._calculateconstants()
        
    def set_stress(self, stress):
        """Set the applied stress.

        Must be a symmetric 3x3 tensor, a 6-vector representing a symmetric
        3x3 tensor, or a number representing the pressure.
        """
        if type(stress) == type(1.0) or type(stress) == type(1):
            stress = array((-stress, -stress, -stress, 0.0, 0.0, 0.0))
        elif stress.shape == (3,3):
            if not self._issymmetric(stress):
                raise ValueError, "The external stress must be a symmetric tensor."
            stress = array((stress[0,0], stress[1,1], stress[2,2], stress[1,2],
                            stress[0,2], stress[0,1]))
        elif stress.shape != (6,):
            raise ValueError, "The external stress has the wrong shape."
        self.externalstress = stress

    def set_mask(self, mask):
        """Set the mask indicating dynamic elements of the computational box.

        If set to None, all elements may change.  If set to a 3-vector
        of ones and zeros, elements which are zero specify directions
        along which the size of the computational box cannot change.
        For example, if mask = {1,1,0} the length of the system along
        the z-axis cannot change, although xz and yz shear is still
        possible.  To disable shear globally, set the mode to diagonal
        (not yet implemented).
        """
        if mask is None:
            mask = ones((3,))
        if not hasattr(mask, "shape"):
            mask = array(mask)        
        if mask.shape != (3,) and mask.shape != (3,3):
            raise "The mask has the wrong shape (must be a 3-vector or 3x3 matrix)"
        else:
            mask = not_equal(mask, 0)  # Make sure it is 0/1

        if mask.shape == (3,):
            self.mask = outer(mask, mask)
        else:
            self.mask = mask
        
    def set_fraction_traceless(self, fracTraceless):
        """set what fraction of the traceless part of the force
        on eta is kept.

        By setting this to zero, the volume may change but the shape may not.
        """
        self.frac_traceless = fracTraceless

    def get_strain_rate(self):
        "Get the strain rate as an upper-triangular 3x3 matrix"
        return array(self.eta, copy=1)

    def set_strain_rate(self, rate):
        "Set the strain rate.  Must be an upper triangular 3x3 matrix."
        if not (rate.shape == (3,3) and self._isuppertriangular(rate)):
            raise ValueError, "Strain rate must be an upper triangular matrix."
        self.eta = rate
        if self.initialized:
            # Recalculate h_past and eta_past so they match the current value.
            self._initialize_eta_h()

    def get_time(self):
        "Get the elapsed time."
        return self.timeelapsed
    
    def run(self, steps):
        """Perform a number of time steps."""
        if not self.initialized:
            self.initialize()
        else:
            if self.have_the_atoms_been_changed():
                raise NotImplementedError, "You have modified the atoms since the last timestep."

        for i in xrange(steps):
            self.step()
            self.nsteps += 1
            self.call_observers()

    def have_the_atoms_been_changed(self):
        "Checks if the user has modified the positions or momenta of the atoms"
        limit = 1e-10
        h = self._getbox()
        if max(abs((h - self.h).ravel())) > limit:
            self._warning("The computational box has been modified.")
            return 1
        expected_r = dot(self.q + 0.5, h)
        err = max(abs((expected_r - self.atoms.get_positions()).ravel())) 
        if err > limit:
            self._warning("The atomic positions have been modified: "+ str(err))
            return 1
        return 0
    
    def step(self):
        """Perform a single time step.
        
        Assumes that the forces and stresses are up to date, and that
        the positions and momenta have not been changed since last
        timestep.
        """
        
        ## Assumes the following variables are OK
        # q_past, q, q_future, p, eta, eta_past, zeta, zeta_past, h, h_past
        #
        # q corresponds to the current positions
        # p must be equal to self.atoms.GetCartesianMomenta()
        # h must be equal to self.atoms.GetUnitCell()
        #
        #print "Making a timestep"
        dt = self.dt
        h_future = self.h_past + 2*dt * dot(self.h, self.eta)
        if self.pfactor_given is None:
            deltaeta = zeros(6, float)
        else:
            stress = self.stresscalculator()
            deltaeta = -2*dt * (self.pfact * linalg.det(self.h)
                                * (stress - self.externalstress))
        
        if self.frac_traceless == 1:
            eta_future = self.eta_past + self.mask * self._makeuppertriangular(deltaeta)
        else:
            trace_part, traceless_part = self._separatetrace(self._makeuppertriangular(deltaeta))
            eta_future = self.eta_past + trace_part + self.frac_traceless * traceless_part

        deltazeta = 2*dt*self.tfact * (self.atoms.get_kinetic_energy()
                                       - self.desiredEkin)
        zeta_future = self.zeta_past + deltazeta
        # Advance time
        #print "Max change in scaled positions:", max(abs(self.q_future.flat - self.q.flat))
        #print "Max change in basis set", max(abs((h_future - self.h).flat))
        self.timeelapsed += dt
        self.h_past = self.h
        self.h = h_future
        self.inv_h = linalg.inv(self.h)
        # Do not throw away the q arrays, they are "magical" on parallel
        # simulations (the contents migrate along with the atoms).
        (self.q_past, self.q, self.q_future) = (self.q, self.q_future,
                                                self.q_past)
        self._setbox_and_positions(self.h,self.q)
        self.eta_past = self.eta
        self.eta = eta_future
        self.zeta_past = self.zeta
        self.zeta = zeta_future
        self._synchronize()  # for parallel simulations.
        self.zeta_integrated += dt * self.zeta
        force = self.forcecalculator()
        # The periodic boundary conditions may have moved the atoms.
        self.post_pbc_fix(fixfuture=0)  
        self._calculate_q_future(force)
        self.atoms.set_momenta(dot(self.q_future-self.q_past, self.h/(2*dt)) *
                               self._getmasses())
        #self.stresscalculator()
        
    def forcecalculator(self):
        return self.atoms.get_forces()
    
    def stresscalculator(self):
        return self.atoms.get_stress()

    def initialize(self):
        """Initialize the dynamics.

        The dynamics requires positions etc for the two last times to
        do a timestep, so the algorithm is not self-starting.  This
        method performs a 'backwards' timestep to generate a
        configuration before the current.
        """
        #print "Initializing the NPT dynamics."
        dt = self.dt
        atoms = self.atoms
        self.h = self._getbox()
        if not self._isuppertriangular(self.h):
            print "I am", self
            print "self.h:"
            print self.h
            print "Min:", min((self.h[1,0], self.h[2,0], self.h[2,1]))
            print "Max:", max((self.h[1,0], self.h[2,0], self.h[2,1]))
            raise NotImplementedError, "Can (so far) only operate on lists of atoms where the computational box is an upper triangular matrix."
        self.inv_h = linalg.inv(self.h)
        # The contents of the q arrays should migrate in parallel simulations.
        self._make_special_q_arrays()
        self.q[:] = dot(self.atoms.get_positions(),
                                self.inv_h) - 0.5
        # zeta and eta were set in __init__
        self._initialize_eta_h()
        deltazeta = dt * self.tfact * (atoms.get_kinetic_energy() -
                                       self.desiredEkin)
        self.zeta_past = self.zeta - deltazeta
        self._calculate_q_past_and_future()
        self.initialized = 1

    def get_gibbs_free_energy(self):
        """Return the Gibb's free energy, which is supposed to be conserved.

        Requires that the energies of the atoms are up to date.

        This is mainly intended as a diagnostic tool.  If called before the
        first timestep, Initialize will be called.
        """
        if not self.initialized:
            self.initialize()
        n = self._getnatoms()
        #tretaTeta = sum(diagonal(matrixmultiply(transpose(self.eta),
        #                                        self.eta)))
        contractedeta = sum((self.eta*self.eta).ravel())
        gibbs = (self.atoms.get_potential_energy() +
                 self.atoms.get_kinetic_energy()
                 - sum(self.externalstress[0:3]) * linalg.det(self.h) / 3.0)
        if self.ttime is not None:
            gibbs += (1.5 * n * self.temperature * (self.ttime * self.zeta)**2
                      + 3 * self.temperature * (n-1) * self.zeta_integrated)
        else:
            assert self.zeta == 0.0
        if self.pfactor_given is not None:
            gibbs += 0.5 / self.pfact * contractedeta
        else:
            assert contractedeta == 0.0
        return gibbs

    def get_center_of_mass_momentum(self):
        "Get the center of mass momentum."
        return self.atoms.get_momenta().sum(0)

    def zero_center_of_mass_momentum(self, verbose=0):
        "Set the center of mass momentum to zero."
        cm = self.get_center_of_mass_momentum()
        abscm = sqrt(sum(cm*cm))
        if verbose and abscm > 1e-4:
            self._warning(self.classname+": Setting the center-of-mass momentum to zero (was %.6g %.6g %.6g)" % tuple(cm))
        self.atoms.set_momenta(self.atoms.get_momenta()
                               - cm / self._getnatoms())
    
    def post_pbc_fix(self, fixfuture=1):
        """Correct for atoms moved by the boundary conditions.

        If the fixfuture argument is 1 (the default), q_future is also
        corrected.  This is not necessary when post_pbc_fix() is called from
        within Timestep(), but must be done when the user calls post_pbc_fix
        (for example if a CNA calculation may have triggered a migration).
        """
        q = dot(self.atoms.get_positions(),
                           self.inv_h) - 0.5
        delta_q = floor(0.5 + (q - self.q))
        self.q += delta_q
        self.q_past += delta_q
        if fixfuture:
            self.q_future += delta_q
        
    def attach_atoms(self, atoms):
        """Assign atoms to a restored dynamics object.

        This function must be called to set the atoms immediately after the
        dynamics object has been read from a trajectory.
        """
        try:
            self.atoms
        except AttributeError:
            pass
        else:
            raise RuntimeError, "Cannot call attach_atoms on a dynamics which already has atoms."
        MolecularDynamics.__init__(self, atoms, self.dt)
        ####self.atoms = atoms
        limit = 1e-6
        h = self._getbox()
        if max(abs((h - self.h).ravel())) > limit:
            raise RuntimeError, "The unit cell of the atoms does not match the unit cell stored in the file."
        self.inv_h = linalg.inv(self.h)
        self._make_special_q_arrays()
        self.q[:] = dot(self.atoms.get_positions(),
                                           self.inv_h) - 0.5
        self._calculate_q_past_and_future()
        self.initialized = 1
        
    def attach(self, function, interval=1, *args, **kwargs):
        """Attach callback function or trajectory.

        At every *interval* steps, call *function* with arguments
        *args* and keyword arguments *kwargs*.
        
        If *function* is a trajectory object, its write() method is
        attached, but if *function* is a BundleTrajectory (or another
        trajectory supporting set_extra_data(), said method is first
        used to instruct the trajectory to also save internal
        data from the NPT dynamics object.
        """
        if hasattr(function, "set_extra_data"):
            # We are attaching a BundleTrajectory or similar
            function.set_extra_data("npt_init",
                                    WeakMethodWrapper(self, "get_init_data"),
                                    once=True)
            function.set_extra_data("npt_dynamics",
                                    WeakMethodWrapper(self, "get_data"))
        MolecularDynamics.attach(self, function, interval, *args, **kwargs)

    def get_init_data(self):
        "Return the data needed to initialize a new NPT dynamics."
        return {'dt': self.dt,
                'temperature': self.temperature,
                'desiredEkin': self.desiredEkin,
                'externalstress': self.externalstress,
                'mask': self.mask,
                'ttime': self.ttime,
                'tfact': self.tfact,
                'pfactor_given': self.pfactor_given,
                'pfact': self.pfact,
                'frac_traceless': self.frac_traceless}
        
    def get_data(self):
        "Return data needed to restore the state."
        return {'eta': self.eta,
                'eta_past': self.eta_past,
                'zeta': self.zeta,
                'zeta_past': self.zeta_past,
                'zeta_integrated': self.zeta_integrated,
                'h': self.h,
                'h_past': self.h_past,
                'timeelapsed': self.timeelapsed}
        
    @classmethod
    def read_from_trajectory(cls, trajectory, frame=-1, atoms=None):
        """Read dynamics and atoms from trajectory (Class method).
        
        Simultaneously reads the atoms and the dynamics from a BundleTrajectory,
        including the internal data of the NPT dynamics object (automatically
        saved when attaching a BundleTrajectory to an NPT object).
        
        Arguments::
        
        trajectory 
            The filename or an open BundleTrajectory object.
        
        frame (optional)
            Which frame to read.  Default: the last.
            
        atoms (optional, internal use only)
            Pre-read atoms.  Do not use. 
        """
        if isinstance(trajectory, str):
            if trajectory.endswith('/'):
                trajectory = trajectory[:-1]
            if trajectory.endswith('.bundle'):
                from ase.io.bundletrajectory import BundleTrajectory
                trajectory = BundleTrajectory(trajectory)
            else:
                raise ValueError("Cannot open '%': unsupported file format" % trajectory)
        # trajectory is now a BundleTrajectory object (or compatible)
        if atoms is None:
            atoms = trajectory[frame]
        init_data = trajectory.read_extra_data('npt_init', 0)
        frame_data = trajectory.read_extra_data('npt_dynamics', frame)
        dyn = cls(atoms, timestep=init_data['dt'], 
                  temperature=init_data['temperature'],
                  externalstress=init_data['externalstress'],
                  ttime=init_data['ttime'],
                  pfactor=init_data['pfactor_given'],
                  mask=init_data['mask'])
        dyn.desiredEkin = init_data['desiredEkin']
        dyn.tfact = init_data['tfact']
        dyn.pfact = init_data['pfact']
        dyn.frac_traceless = init_data['frac_traceless']
        for k, v in frame_data.items():
            setattr(dyn, k, v)
        return (dyn, atoms)
        
    def _getbox(self):
        "Get the computational box."
        return self.atoms.get_cell()

    def _getmasses(self):
        "Get the masses as an Nx1 array."
        return reshape(self.atoms.get_masses(), (-1,1))
    
#    def _getcartesianpositions(self):
#        "Get the cartesian positions of the atoms"
#        return self.atoms.get_positions()
    
#    def _getmomenta(self):
#        "Get the (cartesian) momenta of the atoms"
#        return self.atoms.GetCartesianMomenta()

#    def _getforces(self):
#        "Get the (cartesian) forces of the atoms"
#        return self.atoms.GetCartesianForces()

#    def _setmomenta(self, momenta):
#        "Set the (cartesian) momenta of the atoms"
#        self.atoms.SetCartesianMomenta(momenta)
        
    def _separatetrace(self, mat):
        """return two matrices, one proportional to the identity
        the other traceless, which sum to the given matrix
        """
        tracePart = ((mat[0][0] + mat[1][1] + mat[2][2]) / 3.) * identity(3)
        return tracePart, mat - tracePart

    # A number of convenient helper methods
    def _warning(self, text):
        "Emit a warning."
        sys.stderr.write("WARNING: "+text+"\n")
        sys.stderr.flush()
    
    def _calculate_q_future(self, force):
        "Calculate future q.  Needed in Timestep and Initialization."
        dt = self.dt
        id3 = identity(3)
        alpha = (dt * dt) * dot(force / self._getmasses(),
                                self.inv_h)
        beta = dt * dot(self.h, dot(self.eta + 0.5 * self.zeta * id3,
                                    self.inv_h))
        inv_b = linalg.inv(beta + id3)
        self.q_future[:] = dot(2*self.q + dot(self.q_past, beta - id3) + alpha,
                               inv_b)

    def _calculate_q_past_and_future(self):
        def ekin(p, m = self.atoms.get_masses()):
            p2 = sum(p*p, -1)
            return 0.5 * sum(p2 / m) / len(m)
        p0 = self.atoms.get_momenta()
        m = self._getmasses()
        e0 = ekin(p0)
        p = array(p0, copy=1)
        dt = self.dt
        for i in range(2):
            self.q_past[:] = self.q - dt * dot(p / m, self.inv_h)
            self._calculate_q_future(self.atoms.get_forces())
            p = dot(self.q_future - self.q_past, self.h/(2*dt)) * m
            e = ekin(p)
            if e < 1e-5:
                # The kinetic energy and momenta are virtually zero
                return
            p = (p0 - p) + p0

    def _initialize_eta_h(self):
        self.h_past = self.h - self.dt * dot(self.h, self.eta)
        if self.pfactor_given is None:
            deltaeta = zeros(6, float)
        else:
            deltaeta = (-self.dt * self.pfact * linalg.det(self.h)
                        * (self.atoms.get_stress() - self.externalstress))
        if self.frac_traceless == 1:
            self.eta_past = self.eta - self.mask * self._makeuppertriangular(deltaeta)
        else:
            trace_part, traceless_part = self._separatetrace(self._makeuppertriangular(deltaeta))
            self.eta_past = self.eta - trace_part - self.frac_traceless * traceless_part
        
    
    def _makeuppertriangular(self, sixvector):
        "Make an upper triangular matrix from a 6-vector."
        return array(((sixvector[0], sixvector[5], sixvector[4]),
                      (0,            sixvector[1], sixvector[3]),
                      (0,            0,            sixvector[2])))

    def _isuppertriangular(self, m):
        "Check that a matrix is on upper triangular form."
        return m[1,0] == m[2,0] == m[2,1] == 0.0
    
    def _calculateconstants(self):
        "(Re)calculate some constants when pfactor, ttime or temperature have been changed."
        n = self._getnatoms()
        if self.ttime is None:
            self.tfact = 0.0
        else:
            self.tfact = 2.0 / (3 * n * self.temperature *
                                self.ttime * self.ttime)
        if self.pfactor_given is None:
            self.pfact = 0.0
        else:
            self.pfact = 1.0 / (self.pfactor_given
                                * linalg.det(self._getbox()))
            #self.pfact = 1.0/(n * self.temperature * self.ptime * self.ptime)
        self.desiredEkin = 1.5 * (n - 1) * self.temperature

    def _setbox_and_positions(self, h, q):
        """Set the computational box and the positions."""
        self.atoms.set_cell(h, scale_atoms=True)
        r = dot(q + 0.5, h)
        self.atoms.set_positions(r)

    # A few helper methods, which have been placed in separate methods
    # so they can be replaced in the parallel version.
    def _synchronize(self):
        """Synchronizes eta, h and zeta on all processors in a parallel simulation.

        In a parallel simulation, eta, h and zeta are communicated
        from the master to all slaves, to prevent numerical noise from
        causing them to diverge.

        In a serial simulation, do nothing.
        """
        pass  # This is a serial simulation object.  Do nothing.
    
    def _getnatoms(self):
        """Get the number of atoms.

        In a parallel simulation, this is the total number of atoms on all
        processors.
        """
        return len(self.atoms)
    
    def _make_special_q_arrays(self):
        """Make the arrays used to store data about the atoms.

        In a parallel simulation, these are migrating arrays.  In a
        serial simulation they are ordinary Numeric arrays.
        """
        natoms = len(self.atoms)
        self.q = zeros((natoms,3), float)
        self.q_past = zeros((natoms,3), float)
        self.q_future = zeros((natoms,3), float)

class WeakMethodWrapper:
    """A weak reference to a method.
    
    Create an object storing a weak reference to an instance and 
    the name of the method to call.  When called, calls the method.
    
    Just storing a weak reference to a bound method would not work,
    as the bound method object would go away immediately.
    """
    def __init__(self, obj, method):
        self.obj = weakref.proxy(obj)
        self.method = method
        
    def __call__(self, *args, **kwargs):
        m = getattr(self.obj, self.method)
        return m(*args, **kwargs)

# class _HooverNPTTrajectory:
#     """A Trajectory-like object storing data in a HooverNPT object."""
#     def InitForWrite(self):
#         """Does initialization related to write mode."""
#         self.CreateDimension('unlim', None)
#         self.nc.history = 'ASE NPT trajectory'
#         self.nc.version = '0.1'
#         self.nc.classname = self.atoms.classname
#         self.unlim = 0
#         self.nc.lengthunit = units.GetLengthUnit()
#         self.nc.energyunit = units.GetEnergyUnit()
#         self.conversion = (1, 1)

#     def InitForWriteOrAppend(self):
#         """Does initialization related to write and append mode.

#         Either InitForWrite or InitForReadOrAppend will have been
#         called before calling this method.
#         """
#         names = copy.copy(self.known_names)
#         if self.atoms.ttime is None:
#             del names['ttime']
#         if self.atoms.pfactor_given is None:
#             del names['pfactor_given']
#         for d in names.keys():
#             def getdata(atoms=self.atoms, name=d):
#                 return getattr(atoms, name)
#             self.Add(d, data = getdata)
                     
#     known_names = {
#         #    name                 shape        typecode  once    units
#         # ----------------------------------------------------------------
#         'dt':              ((),                Float,    True,   (1, -0.5)),
#         'temperature':     ((),                Float,    True,   (0, 1)),
#         'desiredEkin':     ((),                Float,    True,   (0, 1)),
#         'externalstress':  ((6,),              Float,    True,   (-3, 1)),
#         'mask':            ((3, 3),            Float,    True,   (0, 0)),
#         'ttime':           ((),                Float,    True,   (1, -0.5)),
#         'tfact':           ((),                Float,    True,   (-2, 0)),
#         'pfactor_given':   ((),                Float,    True,   (-1, 0)),
#         'pfact':           ((),                Float,    True,   (-2, 0)),
#         'frac_traceless':  ((),                Float,    True,   (0, 0)),
#         'eta':             ((3, 3),            Float,    False,  (-1, 0.5)),
#         'eta_past':        ((3, 3),            Float,    False,  (-1, 0.5)),
#         'zeta':            ((),                Float,    False,  (-1, 0.5)),
#         'zeta_past':       ((),                Float,    False,  (-1, 0.5)),
#         'zeta_integrated': ((),                Float,    False,  (0, 0)),
#         'h':               ((3, 3),            Float,    False,  (1, 0)),
#         'h_past':          ((3, 3),            Float,    False,  (1, 0)),
#         'timeelapsed':     ((),                Float,    False,  (1, -0.5))
#         }

#     # This trajectory does not store a list of atoms
#     def GetListOfAtoms(self, frame=None):
#         raise AttributeError, "GetListOfAtoms makes no sense in a HooverNPTTrajectory"

#     # Instead, we store a dynamics
#     def GetDynamics(self, frame=None):
#         """Get a HooverNPT Dynamics object.

#         If a frame number is not given, the current frame is used.

#         The variant of the object (ASE HooverNPT, ASAP Serial/Parallel NPT)
#         will be the same as the stored object.

#         After getting the dynamics, the atoms should be attached with the
#         dynamics.attach_atoms(atoms) method.        
#         """
#         # Bypass calling the normal constructor
#         class Dummy:
#             pass
#         dyn = Dummy()
#         dyn.__class__ = self.getClass(self.nc.classname)
#         vars = self.nc.variables
#         for q in self.known_names.keys():
#             if vars.has_key(q):
#                 once = self.known_names[q][2]
#                 if once:
#                     setattr(dyn, q, vars[q].getValue())
#                 else:
#                     setattr(dyn, q, vars[q][frame])
#         return dyn

#     def getClass(self, classname):
#         "Internal function: turns a class name into a class object."
#         if self.nc.classname == "HooverNPT":
#             return HooverNPT
#         else:
#             raise RuntimeError, ("Cannot create a dynamics of type "
#                                  + self.nc.classname)

# class HooverNPTTrajectory(_HooverNPTTrajectory,NetCDFTrajectory):
#     """A Trajectory-like object storing data in a HooverNPT object."""
#     def __init__(self, filename, dynamics=None, mode=None, interval=1):
#         """Open the NetCDF file.

#         If there is no ``dynamics`` argument, then the file is opened
#         in read mode - otherwise, write or append mode is used.  The
#         ``interval`` argument determines how often the configurations
#         are written to file."""
#         # Call the original constructor, but passing the dynamics instead of
#         # the atoms.
#         if dynamics is not None:
#             # Prevents a circular reference when the trajectory is attached
#             # to the dynamics it observes.
#             dynamics = weakref.proxy(dynamics)
#         NetCDFTrajectory.__init__(self, filename,
#                                   atoms=dynamics,
#                                   mode=mode, interval=interval)

    

