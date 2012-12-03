"""Minimum mode follower for finding saddle points in an unbiased way.

There is, currently, only one implemented method: The Dimer method.

"""

import sys
import time
import warnings

import numpy as np

from ase.optimize.optimize import Optimizer
from math import cos, sin, atan, tan, degrees, pi, sqrt
from ase.parallel import rank, size, world
from ase.calculators.singlepoint import SinglePointCalculator

# Handy vector methods
norm = np.linalg.norm
def normalize(vector):
    """Create a unit vector along *vector*"""
    return vector / norm(vector)
def parallel_vector(vector, base):
    """Extract the components of *vector* that are parallel to *base*"""
    return np.vdot(vector, base) * base
def perpendicular_vector(vector, base):
    """Remove the components of *vector* that are parallel to *base*"""
    return vector - parallel_vector(vector, base)
def rotate_vectors(v1i, v2i, angle):
    """Rotate vectors *v1i* and *v2i* by *angle*"""
    cAng = cos(angle)
    sAng = sin(angle)
    v1o = v1i * cAng + v2i * sAng
    v2o = v2i * cAng - v1i * sAng
    # Ensure the length of the input and output vectors is equal
    return normalize(v1o) * norm(v1i), normalize(v2o) * norm(v2i)

class DimerEigenmodeSearch:
    """An implementation of the Dimer's minimum eigenvalue mode search.

    This class implements the rotational part of the dimer saddle point
    searching method.

    Parameters:

    atoms : MinModeAtoms object
        MinModeAtoms is an extension to the Atoms object, which includes
        information about the lowest eigenvalue mode.
    control : DimerControl object
        Contains the parameters necessary for the eigenmode search.
        If no control object is supplied a default DimerControl
        will be created and used.
    basis : list of xyz-values
        Eigenmode. Must be an ndarray of shape (n, 3).
        It is possible to constrain the eigenmodes to be orthogonal
        to this given eigenmode.

    Notes:

    The code is inspired, with permission, by code written by the Henkelman
    group, which can be found at http://theory.cm.utexas.edu/vtsttools/code/

    References:

    .. [1] Henkelman and Jonsson, JCP 111, 7010 (1999)
    .. [2] Olsen, Kroes, Henkelman, Arnaldsson, and Jonsson, JCP 121,
           9776 (2004).
    .. [3] Heyden, Bell, and Keil, JCP 123, 224101 (2005).
    .. [4] Kastner and Sherwood, JCP 128, 014106 (2008).

    """
    def __init__(self, atoms, control=None, eigenmode=None, basis=None, \
                 **kwargs):
        if hasattr(atoms, 'get_eigenmode'):
            self.atoms = atoms
        else:
            e = 'The atoms object must be a MinModeAtoms object'
            raise TypeError(e)
        self.basis = basis

        if eigenmode is None:
            self.eigenmode = self.atoms.get_eigenmode()
        else:
            self.eigenmode = eigenmode

        if control is None:
            self.control = DimerControl(**kwargs)
            w = 'Missing control object in ' + self.__class__.__name__ + \
                '. Using default: DimerControl()'
            warnings.warn(w, UserWarning)
            if self.control.logfile is not None:
                self.control.logfile.write('DIM:WARN: ' + w + '\n')
                self.control.logfile.flush()
        else:
            self.control = control
            # kwargs must be empty if a control object is supplied
            for key in kwargs:
                e = '__init__() got an unexpected keyword argument \'%s\'' % \
                    (key)
                raise TypeError(e)

        self.dR = self.control.get_parameter('dimer_separation')
        self.logfile = self.control.get_logfile()

    def converge_to_eigenmode(self):
        """Perform an eigenmode search."""
        self.set_up_for_eigenmode_search()
        stoprot = False

        # Load the relevant parameters from control
        f_rot_min = self.control.get_parameter('f_rot_min')
        f_rot_max = self.control.get_parameter('f_rot_max')
        trial_angle = self.control.get_parameter('trial_angle')
        max_num_rot = self.control.get_parameter('max_num_rot')
        extrapolate = self.control.get_parameter('extrapolate_forces')

        while not stoprot:
            if self.forces1E == None:
                self.update_virtual_forces()
            else:
                self.update_virtual_forces(extrapolated_forces = True)
            self.forces1A = self.forces1
            self.update_curvature()
            f_rot_A = self.get_rotational_force()

            # Pre rotation stop criteria
            if norm(f_rot_A) <= f_rot_min:
                self.log(f_rot_A, None)
                stoprot = True
            else:
                n_A = self.eigenmode
                rot_unit_A = normalize(f_rot_A)

                # Get the curvature and its derivative
                c0 = self.get_curvature()
                c0d = np.vdot((self.forces2 - self.forces1), rot_unit_A) / \
                      self.dR

                # Trial rotation (no need to store the curvature)
                # NYI variable trial angles from [3]
                n_B, rot_unit_B = rotate_vectors(n_A, rot_unit_A, trial_angle)
                self.eigenmode = n_B
                self.update_virtual_forces()
                self.forces1B = self.forces1

                # Get the curvature's derivative
                c1d = np.vdot((self.forces2 - self.forces1), rot_unit_B) / \
                      self.dR

                # Calculate the Fourier coefficients
                a1 = c0d * cos(2 * trial_angle) - c1d / \
                     (2 * sin(2 * trial_angle))
                b1 = 0.5 * c0d
                a0 = 2 * (c0 - a1)

                # Estimate the rotational angle
                rotangle = atan(b1 / a1) / 2.0

                # Make sure that you didn't find a maximum
                cmin = a0 / 2.0 + a1 * cos(2 * rotangle) + \
                       b1 * sin(2 * rotangle)
                if c0 < cmin:
                    rotangle += pi / 2.0

                # Rotate into the (hopefully) lowest eigenmode
                # NYI Conjugate gradient rotation
                n_min, dummy = rotate_vectors(n_A, rot_unit_A, rotangle)
                self.update_eigenmode(n_min)

                # Store the curvature estimate instead of the old curvature
                self.update_curvature(cmin)

                self.log(f_rot_A, rotangle)

                # Force extrapolation scheme from [4]
                if extrapolate:
                    self.forces1E = sin(trial_angle - rotangle) / \
                        sin(trial_angle) * self.forces1A + sin(rotangle) / \
                        sin(trial_angle) * self.forces1B + \
                        (1 - cos(rotangle) - sin(rotangle) * \
                        tan(trial_angle / 2.0)) * self.forces0
                else:
                    self.forces1E = None

            # Post rotation stop criteria
            if not stoprot:
                if self.control.get_counter('rotcount') >= max_num_rot:
                    stoprot = True
                elif norm(f_rot_A) <= f_rot_max:
                    stoprot = True

    def log(self, f_rot_A, angle):
        """Log each rotational step."""
        # NYI Log for the trial angle
        if self.logfile is not None:
            if angle:
                l = 'DIM:ROT: %7d %9d %9.4f %9.4f %9.4f\n' % \
                    (self.control.get_counter('optcount'),
                    self.control.get_counter('rotcount'),
                    self.get_curvature(), degrees(angle), norm(f_rot_A))
            else:
                l = 'DIM:ROT: %7d %9d %9.4f %9s %9.4f\n' % \
                    (self.control.get_counter('optcount'),
                    self.control.get_counter('rotcount'),
                    self.get_curvature(), '---------', norm(f_rot_A))
            self.logfile.write(l)
            self.logfile.flush()

    def get_rotational_force(self):
        """Calculate the rotational force that acts on the dimer."""
        rot_force = perpendicular_vector((self.forces1 - self.forces2),
                    self.eigenmode) / (2.0 * self.dR)
        if self.basis is not None:
            if len(self.basis) == len(self.atoms) and len(self.basis[0]) == \
               3 and isinstance(self.basis[0][0], float):
                rot_force = perpendicular_vector(rot_force, self.basis)
            else:
                for base in self.basis:
                    rot_force = perpendicular_vector(rot_force, base)
        return rot_force

    def update_curvature(self, curv = None):
        """Update the curvature in the MinModeAtoms object."""
        if curv:
            self.curvature = curv
        else:
            self.curvature = np.vdot((self.forces2 - self.forces1),
                             self.eigenmode) / (2.0 * self.dR)

    def update_eigenmode(self, eigenmode):
        """Update the eigenmode in the MinModeAtoms object."""
        self.eigenmode = eigenmode
        self.update_virtual_positions()
        self.control.increment_counter('rotcount')

    def get_eigenmode(self):
        """Returns the current eigenmode."""
        return self.eigenmode

    def get_curvature(self):
        """Returns the curvature along the current eigenmode."""
        return self.curvature

    def get_control(self):
        """Return the control object."""
        return self.control

    def update_center_forces(self):
        """Get the forces at the center of the dimer."""
        self.atoms.set_positions(self.pos0)
        self.forces0 = self.atoms.get_forces(real = True)
        self.energy0 = self.atoms.get_potential_energy()

    def update_virtual_forces(self, extrapolated_forces = False):
        """Get the forces at the endpoints of the dimer."""
        self.update_virtual_positions()

        # Estimate / Calculate the forces at pos1
        if extrapolated_forces:
            self.forces1 = self.forces1E.copy()
        else:
            self.forces1 = self.atoms.get_forces(real = True, pos = self.pos1)

        # Estimate / Calculate the forces at pos2
        if self.control.get_parameter('use_central_forces'):
            self.forces2 = 2 * self.forces0 - self.forces1
        else:
            self.forces2 = self.atoms.get_forces(real = True, pos = self.pos2)

    def update_virtual_positions(self):
        """Update the end point positions."""
        self.pos1 = self.pos0 + self.eigenmode * self.dR
        self.pos2 = self.pos0 - self.eigenmode * self.dR

    def set_up_for_eigenmode_search(self):
        """Before eigenmode search, prepare for rotation."""
        self.pos0 = self.atoms.get_positions()
        self.update_center_forces()
        self.update_virtual_positions()
        self.control.reset_counter('rotcount')
        self.forces1E = None

    def set_up_for_optimization_step(self):
        """At the end of rotation, prepare for displacement of the dimer."""
        self.atoms.set_positions(self.pos0)
        self.forces1E = None

class MinModeControl:
    """A parent class for controlling minimum mode saddle point searches.

    Method specific control classes inherit this one. The only thing
    inheriting classes need to implement are the log() method and
    the *parameters* class variable with default values for ALL
    parameters needed by the method in question.
    When instantiating control classes default parameter values can
    be overwritten.

    """
    parameters = {}
    def __init__(self, logfile = '-', eigenmode_logfile=None, **kwargs):
        # Overwrite the defaults with the input parameters given
        for key in kwargs:
            if not key in self.parameters.keys():
                e = 'Invalid parameter >>%s<< with value >>%s<< in %s' % \
                    (key, str(kwargs[key]), self.__class__.__name__)
                raise ValueError(e)
            else:
                self.set_parameter(key, kwargs[key], log = False)

        # Initialize the log files
        self.initialize_logfiles(logfile, eigenmode_logfile)

        # Initialize the counters
        self.counters = {'forcecalls': 0, 'rotcount': 0, 'optcount': 0}

        self.log()

    def initialize_logfiles(self, logfile=None, eigenmode_logfile=None):
        """Set up the log files."""
        # Set up the regular logfile
        if rank != 0:
            logfile = None
        elif isinstance(logfile, str):
            if logfile == '-':
                logfile = sys.stdout
            else:
                logfile = open(logfile, 'a')
        self.logfile = logfile

        # Set up the eigenmode logfile
        if eigenmode_logfile:
            if rank != 0:
                eigenmode_logfile = None
            elif isinstance(eigenmode_logfile, str):
                if eigenmode_logfile == '-':
                    eigenmode_logfile = sys.stdout
                else:
                    eigenmode_logfile = open(eigenmode_logfile, 'a')
        self.eigenmode_logfile = eigenmode_logfile

    def log(self, parameter=None):
        """Log the parameters of the eigenmode search."""
        pass

    def set_parameter(self, parameter, value, log=True):
        """Change a parameter's value."""
        if not parameter in self.parameters.keys():
            e = 'Invalid parameter >>%s<< with value >>%s<<' % \
                (parameter, str(value))
            raise ValueError(e)
        self.parameters[parameter] = value
        if log:
            self.log(parameter)

    def get_parameter(self, parameter):
        """Returns the value of a parameter."""
        if not parameter in self.parameters.keys():
            e = 'Invalid parameter >>%s<<' % \
                (parameter)
            raise ValueError(e)
        return self.parameters[parameter]

    def get_logfile(self):
        """Returns the log file."""
        return self.logfile

    def get_eigenmode_logfile(self):
        """Returns the eigenmode log file."""
        return self.eigenmode_logfile

    def get_counter(self, counter):
        """Returns a given counter."""
        return self.counters[counter]

    def increment_counter(self, counter):
        """Increment a given counter."""
        self.counters[counter] += 1

    def reset_counter(self, counter):
        """Reset a given counter."""
        self.counters[counter] = 0

    def reset_all_counters(self):
        """Reset all counters."""
        for key in self.counters.keys():
            self.counters[key] = 0

class DimerControl(MinModeControl):
    """A class that takes care of the parameters needed for a Dimer search.

    Parameters:

    eigenmode_method: str
        The name of the eigenmode search method.
    f_rot_min: float
        Size of the rotational force under which no rotation will be
        performed.
    f_rot_max: float
        Size of the rotational force under which only one rotation will be
        performed.
    max_num_rot: int
        Maximum number of rotations per optimizer step.
    trial_angle: float
        Trial angle for the finite difference estimate of the rotational
        angle in radians.
    trial_trans_step: float
        Trial step size for the MinModeTranslate optimizer.
    maximum_translation: float
        Maximum step size and forced step size when the curvature is still
        positive for the MinModeTranslate optimizer.
    cg_translation: bool
        Conjugate Gradient for the MinModeTranslate optimizer.
    use_central_forces: bool
        Only calculate the forces at one end of the dimer and extrapolate
        the forces to the other.
    dimer_separation: float
        Separation of the dimer's images.
    initial_eigenmode_method: str
        How to construct the initial eigenmode of the dimer. If an eigenmode
        is given when creating the MinModeAtoms object, this will be ignored.
        Possible choices are: 'gauss' and 'displacement'
    extrapolate_forces: bool
        When more than one rotation is performed, an extrapolation scheme can
        be used to reduce the number of force evaluations.
    displacement_method: str
        How to displace the atoms. Possible choices are 'gauss' and 'vector'.
    gauss_std: float
        The standard deviation of the gauss curve used when doing random
        displacement.
    order: int
        How many lowest eigenmodes will be inverted.
    mask: list of bool
        Which atoms will be moved during displacement.
    displacement_center: int or [float, float, float]
        The center of displacement, nearby atoms will be displaced.
    displacement_radius: float
        When choosing which atoms to displace with the *displacement_center*
        keyword, this decides how many nearby atoms to displace.
    number_of_displacement_atoms: int
        The amount of atoms near *displacement_center* to displace.

    """
    # Default parameters for the Dimer eigenmode search
    parameters = {'eigenmode_method': 'dimer',
                  'f_rot_min': 0.1,
                  'f_rot_max': 1.00,
                  'max_num_rot': 1,
                  'trial_angle': pi / 4.0,
                  'trial_trans_step': 0.001,
                  'maximum_translation': 0.1,
                  'cg_translation': True,
                  'use_central_forces': True,
                  'dimer_separation': 0.0001,
                  'initial_eigenmode_method': 'gauss',
                  'extrapolate_forces': False,
                  'displacement_method': 'gauss',
                  'gauss_std': 0.1,
                  'order': 1,
                  'mask': None, # NB mask should not be a "parameter"
                  'displacement_center': None,
                  'displacement_radius': None,
                  'number_of_displacement_atoms': None}

    # NB: Can maybe put this in EigenmodeSearch and MinModeControl
    def log(self, parameter=None):
        """Log the parameters of the eigenmode search."""
        if self.logfile is not None:
            if parameter is not None:
                l = 'DIM:CONTROL: Updated Parameter: %s = %s\n' % (parameter,
                     str(self.get_parameter(parameter)))
            else:
                l = 'MINMODE:METHOD: Dimer\n'
                l += 'DIM:CONTROL: Search Parameters:\n'
                l += 'DIM:CONTROL: ------------------\n'
                for key in self.parameters:
                    l += 'DIM:CONTROL: %s = %s\n' % (key,
                         str(self.get_parameter(key)))
                l += 'DIM:CONTROL: ------------------\n'
                l += 'DIM:ROT: OPT-STEP ROT-STEP CURVATURE ROT-ANGLE ' + \
                     'ROT-FORCE\n'
            self.logfile.write(l)
            self.logfile.flush()

class MinModeAtoms:
    """Wrapper for Atoms with information related to minimum mode searching.

    Contains an Atoms object and pipes all unknown function calls to that
    object.
    Other information that is stored in this object are the estimate for
    the lowest eigenvalue, *curvature*, and its corresponding eigenmode,
    *eigenmode*. Furthermore, the original configuration of the Atoms
    object is stored for use in multiple minimum mode searches.
    The forces on the system are modified by inverting the component
    along the eigenmode estimate. This eventually brings the system to
    a saddle point.

    Parameters:

    atoms : Atoms object
        A regular Atoms object
    control : MinModeControl object
        Contains the parameters necessary for the eigenmode search.
        If no control object is supplied a default DimerControl
        will be created and used.
    mask: list of bool
        Determines which atoms will be moved when calling displace()
    random_seed: int
        The seed used for the random number generator. Defaults to
        modified version the current time.

    References:

    .. [1] Henkelman and Jonsson, JCP 111, 7010 (1999)
    .. [2] Olsen, Kroes, Henkelman, Arnaldsson, and Jonsson, JCP 121,
           9776 (2004).
    .. [3] Heyden, Bell, and Keil, JCP 123, 224101 (2005).
    .. [4] Kastner and Sherwood, JCP 128, 014106 (2008).

    """
    def __init__(self, atoms, control=None, eigenmodes=None, random_seed=None, **kwargs):
        self.minmode_init = True
        self.atoms = atoms

        # Initialize to None to avoid strange behaviour due to __getattr__
        self.eigenmodes = eigenmodes
        self.curvatures = None

        if control is None:
            self.control = DimerControl(**kwargs)
            w = 'Missing control object in ' + self.__class__.__name__ + \
                '. Using default: DimerControl()'
            warnings.warn(w, UserWarning)
            if self.control.logfile is not None:
                self.control.logfile.write('DIM:WARN: ' + w + '\n')
                self.control.logfile.flush()
        else:
            self.control = control
            logfile = self.control.get_logfile()
            mlogfile = self.control.get_eigenmode_logfile()
            for key in kwargs:
                if key == 'logfile':
                    logfile = kwargs[key]
                elif key == 'eigenmode_logfile':
                    mlogfile = kwargs[key]
                else:
                    self.control.set_parameter(key, kwargs[key])
            self.control.initialize_logfiles(logfile = logfile,
                                             eigenmode_logfile = mlogfile)

        # Seed the randomness
        if random_seed is None:
            t = time.time()
            if size > 1:
                t = world.sum(t) / float(size)
            # Harvest the latter part of the current time
            random_seed = int(('%30.9f' % t)[-9:])
        self.random_state = np.random.RandomState(random_seed)

        # Check the order
        self.order = self.control.get_parameter('order')

        # Construct the curvatures list
        self.curvatures = [100.0] * self.order

        # Save the original state of the atoms.
        self.atoms0 = self.atoms.copy()
        self.save_original_forces()

        # Get a reference to the log files
        self.logfile = self.control.get_logfile()
        self.mlogfile = self.control.get_eigenmode_logfile()

    def save_original_forces(self, force_calculation=False):
        """Store the forces (and energy) of the original state."""
        # NB: Would be nice if atoms.copy() took care of this.
        if self.calc is not None:
            # Hack because some calculators do not have calculation_required
            if (hasattr(self.calc, 'calculation_required') \
               and not self.calc.calculation_required(self.atoms,
               ['energy', 'forces'])) or force_calculation:
                calc = SinglePointCalculator( \
                       self.atoms.get_potential_energy(), \
                       self.atoms.get_forces(), None, None, self.atoms0)
                self.atoms0.set_calculator(calc)

    def initialize_eigenmodes(self, method=None, eigenmodes=None, \
                              gauss_std=None):
        """Make an initial guess for the eigenmode."""
        if eigenmodes is None:
            pos = self.get_positions()
            old_pos = self.get_original_positions()
            if method == None:
                method = \
                     self.control.get_parameter('initial_eigenmode_method')
            if method.lower() == 'displacement' and (pos - old_pos).any():
                eigenmode = normalize(pos - old_pos)
            elif method.lower() == 'gauss':
                self.displace(log = False, gauss_std = gauss_std,
                              method = method)
                new_pos = self.get_positions()
                eigenmode = normalize(new_pos - pos)
                self.set_positions(pos)
            else:
                e = 'initial_eigenmode must use either \'gauss\' or ' + \
                    '\'displacement\', if the latter is used the atoms ' + \
                    'must have moved away from the original positions.' + \
                    'You have requested \'%s\'.' % method
                raise NotImplementedError(e) # NYI
            eigenmodes = [eigenmode]

        # Create random higher order mode guesses
        if self.order > 1:
            if len(eigenmodes) == 1:
                for k in range(1, self.order):
                    pos = self.get_positions()
                    self.displace(log = False, gauss_std = gauss_std,
                                  method = method)
                    new_pos = self.get_positions()
                    eigenmode = normalize(new_pos - pos)
                    self.set_positions(pos)
                    eigenmodes += [eigenmode]

        self.eigenmodes = eigenmodes
        # Ensure that the higher order mode guesses are all orthogonal
        if self.order > 1:
            for k in range(self.order):
                self.ensure_eigenmode_orthogonality(k)
        self.eigenmode_log()

    # NB maybe this name might be confusing in context to
    # calc.calculation_required()
    def calculation_required(self):
        """Check if a calculation is required."""
        return self.minmode_init or self.check_atoms != self.atoms

    def calculate_real_forces_and_energies(self, **kwargs):
        """Calculate and store the potential energy and forces."""
        if self.minmode_init:
            self.minmode_init = False
            self.initialize_eigenmodes(eigenmodes = self.eigenmodes)
        self.rotation_required = True
        self.forces0 = self.atoms.get_forces(**kwargs)
        self.energy0 = self.atoms.get_potential_energy()
        self.control.increment_counter('forcecalls')
        self.check_atoms = self.atoms.copy()

    def get_potential_energy(self):
        """Return the potential energy."""
        if self.calculation_required():
            self.calculate_real_forces_and_energies()
        return self.energy0

    def get_forces(self, real=False, pos=None, **kwargs):
        """Return the forces, projected or real."""
        if self.calculation_required() and pos is None:
            self.calculate_real_forces_and_energies(**kwargs)
        if real and pos is None:
            return self.forces0
        elif real and pos != None:
            old_pos = self.atoms.get_positions()
            self.atoms.set_positions(pos)
            forces = self.atoms.get_forces()
            self.control.increment_counter('forcecalls')
            self.atoms.set_positions(old_pos)
            return forces
        else:
            if self.rotation_required:
                self.find_eigenmodes(order = self.order)
                self.eigenmode_log()
                self.rotation_required = False
                self.control.increment_counter('optcount')
            return self.get_projected_forces()

    def ensure_eigenmode_orthogonality(self, order):
        mode = self.eigenmodes[order - 1].copy()
        for k in range(order - 1):
            mode = perpendicular_vector(mode, self.eigenmodes[k])
        self.eigenmodes[order - 1] = normalize(mode)

    def find_eigenmodes(self, order=1):
        """Launch eigenmode searches."""
        if self.control.get_parameter('eigenmode_method').lower() != 'dimer':
            e = 'Only the Dimer control object has been implemented.'
            raise NotImplementedError(e) # NYI
        for k in range(order):
            if k > 0:
                self.ensure_eigenmode_orthogonality(k + 1)
            search = DimerEigenmodeSearch(self, self.control, \
                eigenmode = self.eigenmodes[k], basis = self.eigenmodes[:k])
            search.converge_to_eigenmode()
            search.set_up_for_optimization_step()
            self.eigenmodes[k] = search.get_eigenmode()
            self.curvatures[k] = search.get_curvature()

    def get_projected_forces(self, pos=None):
        """Return the projected forces."""
        if pos is not None:
            forces = self.get_forces(real = True, pos = pos).copy()
        else:
            forces = self.forces0.copy()

        # Loop through all the eigenmodes
        # NB: Can this be done with a linear combination, instead?
        for k, mode in enumerate(self.eigenmodes):
            # NYI This If statement needs to be overridable in the control
            if self.get_curvature(order = k) > 0.0 and self.order == 1:
                forces = -parallel_vector(forces, mode)
            else:
                forces -= 2 * parallel_vector(forces, mode)
        return forces

    def restore_original_positions(self):
        """Restore the MinModeAtoms object positions to the original state."""
        self.atoms.set_positions(self.get_original_positions())

    def get_barrier_energy(self):
        """The energy difference between the current and original states"""
        try:
            original_energy = self.get_original_potential_energy()
            dimer_energy = self.get_potential_energy()
            return dimer_energy - original_energy
        except RuntimeError:
            w = 'The potential energy is not available, without further ' + \
                'calculations, most likely at the original state.'
            warnings.warn(w, UserWarning)
            return np.nan

    def get_control(self):
        """Return the control object."""
        return self.control

    def get_curvature(self, order='max'):
        """Return the eigenvalue estimate."""
        if order == 'max':
            return max(self.curvatures)
        else:
            return self.curvatures[order - 1]

    def get_eigenmode(self, order=1):
        """Return the current eigenmode guess."""
        return self.eigenmodes[order - 1]

    def get_atoms(self):
        """Return the unextended Atoms object."""
        return self.atoms

    def set_atoms(self, atoms):
        """Set a new Atoms object"""
        self.atoms = atoms

    def set_eigenmode(self, eigenmode, order=1):
        """Set the eigenmode guess."""
        self.eigenmodes[order - 1] = eigenmode

    def set_curvature(self, curvature, order=1):
        """Set the eigenvalue estimate."""
        self.curvatures[order - 1] = curvature

    # Pipe all the stuff from Atoms that is not overwritten.
    # Pipe all requests for get_original_* to self.atoms0.
    def __getattr__(self, attr):
        """Return any value of the Atoms object"""
        if 'original' in attr.split('_'):
            attr = attr.replace('_original_', '_')
            return getattr(self.atoms0, attr)
        else:
            return getattr(self.atoms, attr)

    def displace(self, displacement_vector=None, mask=None, method=None,
                 displacement_center=None, radius=None, number_of_atoms=None,
                 gauss_std=None, mic=True, log=True):
        """Move the atoms away from their current position.

        This is one of the essential parts of minimum mode searches.
        The parameters can all be set in the control object and overwritten
        when this method is run, apart from *displacement_vector*.
        It is preferred to modify the control values rather than those here
        in order for the correct ones to show up in the log file.

        *method* can be either 'gauss' for random displacement or 'vector'
        to perform a predefined displacement.

        *gauss_std* is the standard deviation of the gauss curve that is
        used for random displacement.

        *displacement_center* can be either the number of an atom or a 3D
        position. It must be accompanied by a *radius* (all atoms within it
        will be displaced) or a *number_of_atoms* which decides how many of
        the closest atoms will be displaced.

        *mic* controls the usage of the Minimum Image Convention.

        If both *mask* and *displacement_center* are used, the atoms marked
        as False in the *mask* will not be affected even though they are
        within reach of the *displacement_center*.

        The parameters priority order:
        1) displacement_vector
        2) mask
        3) displacement_center (with radius and/or number_of_atoms)

        If both *radius* and *number_of_atoms* are supplied with
        *displacement_center*, only atoms that fulfill both criteria will
        be displaced.

        """

        # Fetch the default values from the control
        if mask is None:
            mask = self.control.get_parameter('mask')
        if method is None:
            method = self.control.get_parameter('displacement_method')
        if gauss_std is None:
            gauss_std = self.control.get_parameter('gauss_std')
        if displacement_center is None:
            displacement_center = \
                    self.control.get_parameter('displacement_center')
        if radius is None:
            radius = self.control.get_parameter('displacement_radius')
        if number_of_atoms is None:
            number_of_atoms = \
                    self.control.get_parameter('number_of_displacement_atoms')

        # Check for conflicts
        if displacement_vector is not None and method.lower() != 'vector':
            e = 'displacement_vector was supplied but a different method ' + \
                '(\'%s\') was chosen.\n' % str(method)
            raise ValueError(e)
        elif displacement_vector is None and method.lower() == 'vector':
            e = 'A displacement_vector must be supplied when using ' + \
                'method = \'%s\'.\n' % str(method)
            raise ValueError(e)
        elif displacement_center is not None and radius is None and \
           number_of_atoms is None:
            e = 'When displacement_center is chosen, either radius or ' + \
                'number_of_atoms must be supplied.\n'
            raise ValueError(e)

        # Set up the center of displacement mask (c_mask)
        if displacement_center is not None:
            c = displacement_center
            # Construct a distance list
            # The center is an atom
            if type(c) is int:
                # Parse negative indexes
                c = displacement_center % len(self)
                d = [(k, self.get_distance(k, c, mic = mic)) for k in \
                     range(len(self))]
            # The center is a position in 3D space
            elif len(c) == 3 and [type(c_k) for c_k in c] == [float]*3:
                # NB: MIC is not considered.
                d = [(k, norm(self.get_positions()[k] - c)) \
                     for k in range(len(self))]
            else:
                e = 'displacement_center must be either the number of an ' + \
                    'atom in MinModeAtoms object or a 3D position ' + \
                    '(3-tuple of floats).'
                raise ValueError(e)

            # Set up the mask
            if radius is not None:
                r_mask = [dist[1] < radius for dist in d]
            else:
                r_mask = [True for _ in self]

            if number_of_atoms is not None:
                d_sorted = [n[0] for n in sorted(d, key = lambda k: k[1])]
                n_nearest = d_sorted[:number_of_atoms]
                n_mask = [k in n_nearest for k in range(len(self))]
            else:
                n_mask = [True for _ in self]

            # Resolve n_mask / r_mask conflicts
            c_mask = [n_mask[k] and r_mask[k] for k in range(len(self))]
        else:
            c_mask = None

        # Set up a True mask if there is no mask supplied
        if mask is None:
            mask = [True for _ in self]
            if c_mask is None:
                w = 'It was not possible to figure out which atoms to ' + \
                    'displace, Will try to displace all atoms.\n'
                warnings.warn(w, UserWarning)
                if self.logfile is not None:
                    self.logfile.write('MINMODE:WARN: ' + w + '\n')
                    self.logfile.flush()

        # Resolve mask / c_mask conflicts
        if c_mask is not None:
            mask = [mask[k] and c_mask[k] for k in range(len(self))]

        if displacement_vector is None:
            displacement_vector = []
            for k in range(len(self)):
                if mask[k]:
                    diff_line = []
                    for _ in range(3):
                        if method.lower() == 'gauss':
                            if not gauss_std:
                                gauss_std = \
                                self.control.get_parameter('gauss_std')
                            diff = self.random_state.normal(0.0, gauss_std)
                        else:
                            e = 'Invalid displacement method >>%s<<' % \
                                 str(method)
                            raise ValueError(e)
                        diff_line.append(diff)
                    displacement_vector.append(diff_line)
                else:
                    displacement_vector.append([0.0]*3)

        # Remove displacement of masked atoms
        for k in range(len(mask)):
            if not mask[k]:
                displacement_vector[k] = [0.0]*3

        # Perform the displacement and log it
        if log:
            pos0 = self.get_positions()
        self.set_positions(self.get_positions() + displacement_vector)
        if log:
            parameters = {'mask': mask,
                          'displacement_method': method,
                          'gauss_std': gauss_std,
                          'displacement_center': displacement_center,
                          'displacement_radius': radius,
                          'number_of_displacement_atoms': number_of_atoms}
            self.displacement_log(self.get_positions() - pos0, parameters)

    def eigenmode_log(self):
        """Log the eigenmodes (eigenmode estimates)"""
        if self.mlogfile is not None:
            l = 'MINMODE:MODE: Optimization Step: %i\n' % \
                   (self.control.get_counter('optcount'))
            for m_num, mode in enumerate(self.eigenmodes):
                l += 'MINMODE:MODE: Order: %i\n' % m_num
                for k in range(len(mode)):
                    l += 'MINMODE:MODE: %7i %15.8f %15.8f %15.8f\n' % (k,
                         mode[k][0], mode[k][1], mode[k][2])
            self.mlogfile.write(l)
            self.mlogfile.flush()

    def displacement_log(self, displacement_vector, parameters):
        """Log the displacement"""
        if self.logfile is not None:
            lp = 'MINMODE:DISP: Parameters, different from the control:\n'
            mod_para = False
            for key in parameters:
                if parameters[key] != self.control.get_parameter(key):
                    lp += 'MINMODE:DISP: %s = %s\n' % (str(key),
                                                       str(parameters[key]))
                    mod_para = True
            if mod_para:
                l = lp
            else:
                l = ''
            for k in range(len(displacement_vector)):
                l += 'MINMODE:DISP: %7i %15.8f %15.8f %15.8f\n' % (k,
                     displacement_vector[k][0], displacement_vector[k][1],
                     displacement_vector[k][2])
            self.logfile.write(l)
            self.logfile.flush()

    def summarize(self):
        """Summarize the Minimum mode search."""
        if self.logfile is None:
            logfile = sys.stdout
        else:
            logfile = self.logfile

        c = self.control
        label = 'MINMODE:SUMMARY: '

        l = label + '-------------------------\n'
        l += label + 'Barrier: %16.4f\n' % self.get_barrier_energy()
        l += label + 'Curvature: %14.4f\n' % self.get_curvature()
        l += label + 'Optimizer steps: %8i\n' % c.get_counter('optcount')
        l += label + 'Forcecalls: %13i\n' % c.get_counter('forcecalls')
        l += label + '-------------------------\n'

        logfile.write(l)

class MinModeTranslate(Optimizer):
    """An Optimizer specifically tailored to minimum mode following."""
    def __init__(self, atoms, logfile='-', trajectory=None):
        Optimizer.__init__(self, atoms, None, logfile, trajectory)

        self.control = atoms.get_control()

        # Make a header for the log
        if self.logfile is not None:
            l = ''
            if isinstance(self.control, DimerControl):
                l = 'MinModeTranslate: STEP      TIME          ENERGY    ' + \
                    'MAX-FORCE     STEPSIZE    CURVATURE  ROT-STEPS\n'
            self.logfile.write(l)
            self.logfile.flush()

        # Load the relevant parameters from control
        self.cg_on = self.control.get_parameter('cg_translation')
        self.trial_step = self.control.get_parameter('trial_trans_step')
        self.max_step = self.control.get_parameter('maximum_translation')

        # Start conjugate gradient
        if self.cg_on:
            self.cg_init = True

    def initialize(self):
        """Set initial values."""
        self.r0 = None
        self.f0 = None

    def run(self, fmax=0.05, steps=100000000):
        """Run structure optimization algorithm.

        This method will return when the forces on all individual
        atoms are less than *fmax* or when the number of steps exceeds
        *steps*.

        """

        self.fmax = fmax
        step = 0
        while step < steps:
            f = self.atoms.get_forces()
            self.call_observers()
            if self.converged(f):
                self.log(f, None)
                return
            self.step(f)
            self.nsteps += 1
            step += 1

    def step(self, f):
        """Perform the optimization step."""
        atoms = self.atoms
        r = atoms.get_positions()
        curv = atoms.get_curvature()
        f0p = f.copy()
        r0 = r.copy()
        direction = f0p.copy()
        if self.cg_on:
            direction = self.get_cg_direction(direction)
        direction = normalize(direction)
        if curv > 0.0:
            step = direction * self.max_step
        else:
            r0t = r0 + direction * self.trial_step
            f0tp = self.atoms.get_projected_forces(r0t)
            F = np.vdot((f0tp + f0p), direction) / 2.0
            C = np.vdot((f0tp - f0p), direction) / self.trial_step
            step = ( -F / C + self.trial_step / 2.0 ) * direction
            if norm(step) > self.max_step:
                step = direction * self.max_step
        self.log(f0p, norm(step))

        atoms.set_positions(r + step)

        self.f0 = f.flat.copy()
        self.r0 = r.flat.copy()

    def get_cg_direction(self, direction):
        """Apply the Conjugate Gradient algorithm to the step direction."""
        if self.cg_init:
            self.cg_init = False
            self.direction_old = direction.copy()
            self.cg_direction = direction.copy()
        old_norm = np.vdot(self.direction_old, self.direction_old)
        # Polak-Ribiere Conjugate Gradient
        if old_norm != 0.0:
            betaPR = np.vdot(direction, (direction - self.direction_old)) / \
                     old_norm
        else:
            betaPR = 0.0
        if betaPR < 0.0:
            betaPR = 0.0
        self.cg_direction = direction + self.cg_direction * betaPR
        self.direction_old = direction.copy()
        return self.cg_direction.copy()

    def log(self, f, stepsize):
        """Log each step of the optimization."""
        if self.logfile is not None:
            T = time.localtime()
            e = self.atoms.get_potential_energy()
            fmax = sqrt((f**2).sum(axis = 1).max())
            rotsteps = self.atoms.control.get_counter('rotcount')
            curvature = self.atoms.get_curvature()
            l = ''
            if stepsize:
                if isinstance(self.control, DimerControl):
                    l = '%s: %4d  %02d:%02d:%02d %15.6f %12.4f %12.6f ' \
                        '%12.6f %10d\n' % ('MinModeTranslate', self.nsteps,
                         T[3], T[4], T[5], e, fmax, stepsize, curvature,
                         rotsteps)
            else:
                if isinstance(self.control, DimerControl):
                    l = '%s: %4d  %02d:%02d:%02d %15.6f %12.4f %s ' \
                        '%12.6f %10d\n' % ('MinModeTranslate', self.nsteps,
                         T[3], T[4], T[5], e, fmax, '    --------',
                         curvature, rotsteps)
            self.logfile.write(l)
            self.logfile.flush()

def read_eigenmode(mlog, index = -1):
    """Read an eigenmode.
    To access the pre optimization eigenmode set index = 'null'.

    """
    if isinstance(mlog, str):
        f = open(mlog, 'r')
    else:
        f = mlog

    lines = f.readlines()

    # Detect the amount of atoms and iterations
    k = 2
    while lines[k].split()[1].lower() not in ['optimization', 'order']:
        k += 1
    n = k - 2
    n_itr = (len(lines) / (n + 1)) - 2

    # Locate the correct image.
    if type(index) == str:
        if index.lower() == 'null':
            i = 0
        else:
            i = int(index) + 1
    else:
        if index >= 0:
            i = index + 1
        else:
            if index < -n_itr - 1:
                raise IndexError('list index out of range')
            else:
                i = index

    mode = np.ndarray(shape = (n, 3), dtype = float)
    k_atom = 0
    for k in range(1, n + 1):
        line = lines[i * (n + 1) + k].split()
        for k_dim in range(3):
            mode[k_atom][k_dim] = float(line[k_dim + 2])
        k_atom += 1

    return mode

# Aliases
DimerAtoms = MinModeAtoms
DimerTranslate = MinModeTranslate
