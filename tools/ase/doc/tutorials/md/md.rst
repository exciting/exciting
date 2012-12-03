.. _md_tutorial:

==================
Molecular dynamics
==================

.. note::

  These examples *can* be used without Asap installed, then
  the ase.EMT calculator (implemented in Python) is used, but nearly
  superhuman patience is required.

Here we demonstrate now simple molecular dynamics is performed.  A
crystal is set up, the atoms are given momenta corresponding to a
temperature of 300K, then Newtons second law is integrated numerically
with a time step of 5 fs (a good choice for copper).

.. literalinclude:: moldyn1.py

Note how the total energy is conserved, but the kinetic energy quickly
drops to half the expected value.  Why?


Instead of printing within a loop, it is possible to use an "observer"
to observe the atoms and do the printing (or more sophisticated
analysis).

.. literalinclude:: moldyn2.py

Constant temperature MD
=======================

Often, you want to control the temperature of an MD simulation.  This
can be done with the Langevin dynamics module.  In the previous
examples, replace the line `dyn = VelocityVerlet(...)` with::

  dyn = Langevin(atoms, 5*units.fs, T*units.kB, 0.002)

where T is the desired temperature in Kelvin.  You also need to import
Langevin, see the class below.

The Langevin dynamics will then slowly adjust the total energy of the
system so the temperature approaches the desired one.

As a slightly less boring example, let us use this to melt a chunck of
copper by starting the simulation without any momentum of the atoms
(no kinetic energy), and with a desired temperature above the melting
point.  We will also save information about the atoms in a trajectory
file called moldyn3.traj.  

.. literalinclude:: moldyn3.py

After running the simulation, you can study the result with the
command

::

  ag moldyn3.traj

Try plotting the kinetic energy.  You will *not* see a well-defined
melting point due to finite size effects (including surface melting),
but you will probably see an almost flat region where the inside of
the system melts.  The outermost layers melt at a lower temperature.

.. note::

  The Langevin dynamics will by default keep the position and momentum
  of the center of mass unperturbed. This is another improvement over
  just setting momenta corresponding to a temperature, as we did before.


Isolated particle MD
====================

When simulating isolated particles with MD, it is sometimes preferable
to set random momenta corresponding to a spefic temperature and let the
system evolve freely. With a relatively high temperature, the is however
a risk that the collection of atoms will drift out of the simulation box
because the randomized momenta gave the center of mass a small but
non-zero velocity too.

Let us see what happens when we propagate a nanoparticle for a long time:

.. literalinclude:: moldyn4.py

After running the simulation, use :program:`ag` to compare the results
with how it looks if you comment out either the line that says
:epydoc:`ase.md.velocitydistribution.Stationary` (atoms),
:epydoc:`ase.md.velocitydistribution.ZeroRotation` (atoms) or both.

::

  ag moldyn4.traj

Try playing the movie with a high frame rate and set frame skipping to a
low number. Can you spot the subtle difference?
