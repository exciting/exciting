.. module:: ase.io.trajectory
   :synopsis: Trajectory input-output module

================
Trajectory files
================

The :mod:`io.trajectory` module defines Trajectory objects, that is
objects storing the temporal evolution of a simulation.  A Trajectory
file contains one or more :class:`~ase.atoms.Atoms` objects, usually to be interpreted as
a time series, although that is not a requirement.

The :mod:`io.trajectory` module currently defines two kinds of
Trajectory files, the PickleTrajectory and the BundleTrajectory.
PickleTrajectory is the recommended Trajectory format,
BundleTrajectory is only intended for large molecular dynamics
simulations (large meaning millions of atoms).

In the future, other kinds of Trajectories may be defined, with
similar Python interface but with different underlying file formats.

PickleTrajectory
================

The PickleTrajectory has the interface

.. autoclass:: ase.io.trajectory.PickleTrajectory
   :members:

Note that there is apparently no methods for reading the trajectory.
Reading is instead done by indexing the trajectory, or by iterating
over the trajectory: ``traj[0]`` and ``traj[-1]`` return the first and
last :class:`~ase.atoms.Atoms` object in the trajectory.

Examples
--------

Reading a configuration::

    from ase.io.trajectory import PickleTrajectory
    traj = PickleTrajectory("example.traj")
    atoms = traj[-1]

Reading all configurations::

    traj = PickleTrajectory("example.traj")
    for atoms in traj:
        # Analyze atoms

Writing every 100th time step in a molecular dynamics simulation::

    # dyn is the dynamics (e.g. VelocityVerlet, Langevin or similar)
    traj = PickleTrajectory("example.traj", "w", atoms)
    dyn.attach(traj.write, interval=100)
    dyn.run(10000)
    traj.close()


BundleTrajectory
================

The BundleTrajectory has the interface

.. autoclass:: ase.io.bundletrajectory.BundleTrajectory
   :members:


See also
========

The function :func:`ase.io.write` can write a single
:class:`~ase.atoms.Atoms` object to a Trajectory file.

The function :func:`ase.io.read` can read an :class:`~ase.atoms.Atoms`
object from a Trajectory file, per default it reads the last one.



