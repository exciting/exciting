"""Berendsen NVT dynamics class."""

import sys
import numpy as np
from ase.md.md import MolecularDynamics
from ase.parallel import world


class NVTBerendsen(MolecularDynamics):
    """Berendsen (constant N, V, T) molecular dynamics.

    Usage: NVTBerendsen(atoms, timestep, temperature, taut, fixcm)

    atoms
        The list of atoms.
        
    timestep
        The time step.

    temperature
        The desired temperature, in Kelvin.

    taut
        Time constant for Berendsen temperature coupling.

    fixcm
        If True, the position and momentum of the center of mass is
        kept unperturbed.  Default: True.

    """

    def __init__(self, atoms, timestep, temperature, taut, fixcm=True,
                 trajectory=None, logfile=None, loginterval=1,
                 communicator=world):

        MolecularDynamics.__init__(self, atoms, timestep, trajectory, 
                                   logfile, loginterval)
        self.taut = taut
        self.temperature = temperature
        self.fixcm = fixcm  # will the center of mass be held fixed?
        self.communicator = communicator

    def set_taut(self, taut):
        self.taut = taut

    def get_taut(self):
        return self.taut

    def set_temperature(self, temperature):
        self.temperature = temperature

    def get_temperature(self):
        return self.temperature

    def set_timestep(self, timestep):
        self.dt = timestep

    def get_timestep(self):
        return self.dt

    def scale_velocities(self):
        """ Do the NVT Berendsen velocity scaling """
        tautscl = self.dt / self.taut
        old_temperature = self.atoms.get_temperature()

        scl_temperature = np.sqrt(1.0+ (self.temperature/ old_temperature- 1.0)
                                  *tautscl)
        #limit the velocity scaling to reasonable values
        if scl_temperature > 1.1:
            scl_temperature = 1.1
        if scl_temperature < 0.9:
            scl_temperature = 0.9
        
        atoms = self.atoms
        p = self.atoms.get_momenta()
        p = scl_temperature * p 
        self.atoms.set_momenta(p)
        return 


    def step(self, f):
        """ move one timestep forward using Berenden NVT molecular dynamics."""
        self.scale_velocities()

        #one step velocity verlet
        atoms = self.atoms
        p = self.atoms.get_momenta()
        p += 0.5 * self.dt * f

        if self.fixcm:
            # calculate the center of mass
            # momentum and subtract it
            psum = p.sum(axis=0) / float(len(p))
            p = p - psum

        self.atoms.set_positions(self.atoms.get_positions() +
             self.dt * p / self.atoms.get_masses()[:,np.newaxis])

        # We need to store the momenta on the atoms before calculating
        # the forces, as in a parallel Asap calculation atoms may
        # migrate during force calculations, and the momenta need to
        # migrate along with the atoms.  For the same reason, we
        # cannot use self.masses in the line above.

        self.atoms.set_momenta(p)
        f = self.atoms.get_forces()
        atoms.set_momenta(self.atoms.get_momenta() + 0.5 * self.dt * f)

        return f

