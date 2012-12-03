# VelocityDistributions.py -- set up a velocity distribution

"""Module for setting up e.g. Maxwell-Boltzmann velocity distributions.

Currently, only one function is defined, MaxwellBoltzmannDistribution,
which sets the momenta of a list of atoms according to a
Maxwell-Boltzmann distribution at a given temperature.
"""

import sys
import numpy as np
from ase.parallel import world


def _maxwellboltzmanndistribution(masses, temp, communicator=world):
    # For parallel GPAW simulations, the random velocities should be
    # distributed.  Uses gpaw world communicator as default, but allow
    # option of specifying other communicator (for ensemble runs)
    xi = np.random.standard_normal((len(masses), 3))
    communicator.broadcast(xi, 0)
    momenta = xi * np.sqrt(masses * temp)[:, np.newaxis]
    return momenta


def MaxwellBoltzmannDistribution(atoms, temp, communicator=world,
                                 force_temp=False):
    """Sets the momenta to a Maxwell-Boltzmann distribution. temp should be
    fed in energy units; i.e., for 300 K use temp=300.*units.kB. If
    force_temp is set to True, it scales the random momenta such that the
    temperature request is precise.
    """
    momenta = _maxwellboltzmanndistribution(atoms.get_masses(), temp,
                                            communicator)
    atoms.set_momenta(momenta)
    if force_temp:
        temp0 = atoms.get_kinetic_energy() / len(atoms) / 1.5
        gamma = temp / temp0
        atoms.set_momenta(atoms.get_momenta() * np.sqrt(gamma))



def Stationary(atoms):
    "Sets the center-of-mass momentum to zero."
    p = atoms.get_momenta()
    p0 = np.sum(p, 0)
    # We should add a constant velocity, not momentum, to the atoms
    m = atoms.get_masses()
    mtot = np.sum(m)
    v0 = p0 / mtot
    p -= v0 * m[:, np.newaxis]
    atoms.set_momenta(p)


def ZeroRotation(atoms):
    "Sets the total angular momentum to zero by counteracting rigid rotations."
    # Find the principal moments of inertia and principal axes basis vectors   
    Ip, basis = atoms.get_moments_of_inertia(vectors=True)
    # Calculate the total angular momentum and transform to principal basis
    Lp = np.dot(basis, atoms.get_angular_momentum())
    # Calculate the rotation velocity vector in the principal basis, avoiding
    # zero division, and transform it back to the cartesian coordinate system
    omega = np.dot(np.linalg.inv(basis), np.select([Ip > 0], [Lp / Ip]))
    # We subtract a rigid rotation corresponding to this rotation vector
    com = atoms.get_center_of_mass()
    positions = atoms.get_positions()
    positions -= com  # translate center of mass to origin
    velocities = atoms.get_velocities()
    atoms.set_velocities(velocities - np.cross(omega, positions))
