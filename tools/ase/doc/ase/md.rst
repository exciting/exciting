==================
Molecular dynamics
==================

.. module:: md
   :synopsis: Molecular Dynamics

Typical computer simulations involve moving the atoms around, either
to optimize a structure (energy minimization) or to do molecular
dynamics.  This chapter discusses molecular dynamics, energy
minimization algorithms will be discussed in the :mod:`optimize`
section.

A molecular dynamics object will operate on the atoms by moving them
according to their forces - it integrates Newton's second law
numerically.  A typical molecular dynamics simulation will use the
`Velocity Verlet dynamics`_.  You create the
:class:`VelocityVerlet` object, giving it the atoms and a time step, and then
you perform dynamics by calling its :meth:`run` method::

  dyn = VelocityVerlet(atoms, 5.0 * units.fs)
  dyn.run(1000)  # take 1000 steps

A number of different algorithms can be used to perform molecular
dynamics, with slightly different results.  

Choosing the time step
======================

All the dynamics objects need a time step.  Choosing it too small will
waste computer time, choosing it too large will make the dynamics
unstable, typically the energy increases dramatically (the system
"blows up").  If the time step is only a little to large, the lack of
energy conservation is most obvious in `Velocity Verlet dynamics`_,
where energy should otherwise be conserved.

Experience has shown that 5 femtoseconds is a good choice for most metallic
systems.  Systems with light atoms (e.g. hydrogen) and/or with strong
bonds (carbon) will need a smaller time step.

All the dynamics objects documented here are sufficiently related to
have the same optimal time step.


File output
===========

The time evolution of the system can be saved in a trajectory file,
by creating a trajectory object, and attaching it to the dynamics
object.  This is documented in the module :mod:`ase.io.trajectory`. 

Unlike the geometry optimization classes, the molecular dynamics
classes do not support giving a trajectory file name in the
constructor.  Instead the trajectory must be attached explicitly to
the dynamics, and it is *stongly recommended* to use the optional
``interval`` argument, so every time step is not written to the file.


Logging
=======

A logging mechanism is provided, printing time; total, potential and
kinetic energy; and temperature (calculated from the kinetic energy).
It is enabled by giving the ``logfile`` argument when the dynamics
object is created, ``logfile`` may be an open file, a filename or the
string '-' meaning standard output.  Per default, a line is printed
for each timestep, specifying the ``loginterval`` argument will chance
this to a more reasonable frequency.

The logging can be customized by explicitly attaching a
:class:`ase.md.MDLogger` object to the dynamics::

  from ase.md import MDLogger
  dyn = VelocityVerlet(atoms, dt=2*ase.units.fs)
  dyn.attach(MDLogger(dyn, atoms, 'md.log', header=False, stress=False,
             peratom=True, mode="a"), interval=1000)

This example will skip the header line and write energies per atom
instead of total energies.  The parameters are

  ``header``: Print a header line defining the columns.

  ``stress``: Print the six components of the stress tensor.

  ``peratom``:  Print energy per atom instead of total energy.

  ``mode``:  If 'a', append to existing file, if 'w' overwrite
  existing file.

Despite appearances, attaching a logger like this does *not* create a
cyclic reference to the dynamics.

.. note::

   If building your own logging class, be sure not to attach the dynamics
   object directly to the logging object. Instead, create a weak reference
   using the ``proxy`` method of the ``weakref`` package. See the
   `ase.md.MDLogger` source code for an example. (If this is not done, a
   cyclic reference may be created which can cause certain calculators,
   such as Jacapo, to not terminate correctly.)



Constant NVE simulations (the microcanonical ensemble)
======================================================

Newton's second law preserves the total energy of the system, and a
straightforward integration of Newton's second law therefore leads to
simulations preserving the total energy of the system (E), the number
of atoms (N) and the volume of the system (V).  The most appropriate
algorithm for doing this is velocity Verlet dynamics, since it gives
very good long-term stability of the total energy even with quite
large time steps.  Fancier algorithms such as Runge-Kutta may give
very good short-term energy preservation, but at the price of a slow
drift in energy over longer timescales, causing trouble for long
simulations.

In a typical NVE simulation, the temperature will remain approximately
constant, but if significant structural changes occurs they may result
in temperature changes.  If external work is done on the system, the
temperature is likely to rise significantly.

Velocity Verlet dynamics
------------------------

.. module:: md.verlet

.. class:: VelocityVerlet(atoms, timestep)


``VelocityVerlet`` is the only dynamics implementing the NVE ensemble.
It requires two arguments, the atoms and the time step.  Choosing
a too large time step will immediately be obvious, as the energy will
increase with time, often very rapidly.

Example: See the tutorial :ref:`md_tutorial`.



Constant NVT simulations (the canonical ensemble)
=================================================

Since Newton's second law conserves energy and not temperature,
simulations at constant temperature will somehow involve coupling the
system to a heat bath.  This cannot help being somewhat artificial.
Two different approaches are possible within ASE.  In Langevin
dynamics, each atom is coupled to a heat bath through a fluctuating
force and a friction term.  In Nosé-Hoover dynamics, a term
representing the heat bath through a single degree of freedom is
introduced into the Hamiltonian.

Langevin dynamics
-----------------

.. module:: md.langevin

.. class:: Langevin(atoms, timestep, temperature, friction)


The Langevin class implements Langevin dynamics, where a (small)
friction term and a fluctuating force are added to Newton's second law
which is then integrated numerically.  The temperature of the heat
bath and magnitude of the friction is specified by the user, the
amplitude of the fluctuating force is then calculated to give that
temperature.  This procedure has some physical justification: in a
real metal the atoms are (weakly) coupled to the electron gas, and the
electron gas therefore acts like a heat bath for the atoms.  If heat
is produced locally, the atoms locally get a temperature that is
higher than the temperature of the electrons, heat is transferred to
the electrons and then rapidly transported away by them.  A Langevin
equation is probably a reasonable model for this process.

A disadvantage of using Langevin dynamics is that if significant heat
is produced in the simulation, then the temperature will stabilize at
a value higher than the specified temperature of the heat bath, since
a temperature difference between the system and the heat bath is
necessary to get a finite heat flow.  Another disadvantage is that the
fluctuating force is stochastic in nature, so repeating the simulation
will not give exactly the same trajectory.

When the ``Langevin`` object is created, you must specify a time step,
a temperature (in energy units) and a friction.  Typical values for
the friction are 0.01-0.02 atomic units.

::

  # Room temperature simulation
  dyn = Langevin(atoms, 5 * units.fs, units.kB * 300, 0.002)

Both the friction and the temperature can be replaced with arrays
giving per-atom values.  This is mostly useful for the friction, where
one can choose a rather high friction near the boundaries, and set it
to zero in the part of the system where the phenomenon being studied
is located.



Nosé-Hoover dynamics
--------------------

In Nosé-Hoover dynamics, an extra term is added to the Hamiltonian
representing the coupling to the heat bath.  From a pragmatic point of
view one can regard Nosé-Hoover dynamics as adding a friction term to
Newton's second law, but dynamically changing the friction coefficient
to move the system towards the desired temperature.  Typically the
"friction coefficient" will fluctuate around zero.

Nosé-Hoover dynamics is not implemented as a separate class, but is a
special case of NPT dynamics.


Berendsen NVT dynamics
-----------------------
.. module:: md.nvtberendsen

.. class:: NVTBerendsen(atoms, timestep, temperature, taut, fixcm)

In Berendsen NVT simulations the velocities are scaled to achieve the desired 
temperature. The speed of the scaling is determined by the parameter taut.

This method does not result proper NVT sampling but it usually is 
sufficiently good in practise (with large taut). For discussion see 
the gromacs manual at www.gromacs.org.

*atoms*:
    The list of atoms.
    
*timestep*:
    The time step.

*temperature*:
    The desired temperature, in Kelvin.

*taut*:
    Time constant for Berendsen temperature coupling.

*fixcm*:
    If True, the position and momentum of the center of mass is
    kept unperturbed.  Default: True.

::

  # Room temperature simulation (300K, 0.1 fs time step)
  dyn = NVTBerendsen(atoms, 0.1 * units.fs, 300, taut=0.5*1000*units.fs)



Constant NPT simulations (the isothermal-isobaric ensemble)
===========================================================

.. module:: md.npt

.. class:: NPT(atoms, timestep, temperature, externalstress, ttime, pfactor, mask=None) 

Dynamics with constant pressure (or optionally, constant stress) and
constant temperature (NPT or N,stress,T ensemble).  It uses the
combination of Nosé-Hoover and Parrinello-Rahman dynamics proposed by
Melchionna et al. [1] and later modified by Melchionna [2].  The
differential equations are integrated using a centered difference
method [3].  Details of the implementation are available in the
document XXX NPTdynamics.tex, distributed with the module.

The dynamics object is called with the following parameters:

*atoms*:
  The atoms object.

*timestep*:
  The timestep in units matching eV, Å, u.  Use the *units.fs* constant.

*temperature*:
  The desired temperature in eV.

*externalstress*:
  The external stress in eV/Å^3.  Either a symmetric
  3x3 tensor, a 6-vector representing the same, or a scalar
  representing the pressure.  Note that the stress is positive in
  tension whereas the pressure is positive in compression: giving a
  scalar p is equivalent to giving the tensor (-p. -p, -p, 0, 0, 0).

*ttime*:
  Characteristic timescale of the thermostat.  Set to None to
  disable the thermostat.

*pfactor*:
  A constant in the barostat differential equation.  If a
  characteristic barostat timescale of ptime is desired, set pfactor
  to ptime^2 * B (where B is the Bulk Modulus).  Set to None to
  disable the barostat.  Typical metallic bulk moduli are of the order
  of 100 GPa or 0.6 eV/Å^3.

*mask=None*:
  Optional argument.  A tuple of three integers (0 or 1),
  indicating if the system can change size along the three Cartesian
  axes.  Set to (1,1,1) or None to allow a fully flexible
  computational box.  Set to (1,1,0) to disallow elongations along the
  z-axis etc.


Useful parameter values:

* The same *timestep* can be used as in Verlet dynamics, i.e. 5 fs is fine
  for bulk copper.

* The *ttime* and *pfactor* are quite critical[4], too small values may
  cause instabilites and/or wrong fluctuations in T / p.  Too
  large values cause an oscillation which is slow to die.  Good
  values for the characteristic times seem to be 25 fs for *ttime*,
  and 75 fs for *ptime* (used to calculate pfactor), at least for
  bulk copper with 15000-200000 atoms.  But this is not well
  tested, it is IMPORTANT to monitor the temperature and
  stress/pressure fluctuations.

It has the following methods:

.. method:: NPT.run(n):

  Perform n timesteps.

.. method:: NPT.initialize():

  Estimates the dynamic variables for time=-1 to start the
  algorithm.  This is automatically called before the first timestep.

.. method:: NPT.set_stress():

  Set the external stress.  Use with care.  It is
  preferable to set the right value when creating the object.

.. method:: NPT.set_mask():

  Change the mask.  Use with care, as you may "freeze" a
  fluctuation in the strain rate.
  
.. method:: NPT.set_strainrate(eps):

  Set the strain rate.  ``eps`` must be an upper-triangular matrix.
  If you set a strain rate along a direction that is "masked out"
  (see ``set_mask``), the strain rate along that direction will be
  maintained constantly.

.. method:: NPT.get_gibbs_free_energy():

  Gibbs free energy is supposed to be
  preserved by this dynamics.  This is mainly intended as a diagnostic
  tool.

References:

[1] S. Melchionna, G. Ciccotti and B. L. Holian, Molecular Physics
78, p. 533 (1993).

[2] S. Melchionna, Physical Review E 61, p. 6165 (2000).

[3] B. L. Holian, A. J. De Groot, W. G. Hoover, and C. G. Hoover,
Physical Review A 41, p. 4552 (1990).

[4] F. D. Di Tolla and M. Ronchetti, Physical Review E 48, p. 1726 (1993).

.. seealso::
    
   The :term:`API` documentation: :epydoc:`ase.md`


Berendsen NPT dynamics
-----------------------
.. module:: md.nptberendsen

.. class:: NPTBerendsen(atoms, timestep, temperature, taut, fixcm, pressure, taup,compressibility)

In Berendsen NPT simulations the velocities are scaled to achieve the desired 
temperature. The speed of the scaling is determined by the parameter taut.

The atom positions and the simulation cell are scaled in order to achieve 
the desired pressure. 

This method does not result proper NPT sampling but it usually is 
sufficiently good in practise (with large taut and taup). For discussion see 
the gromacs manual at www.gromacs.org. or amber at ambermd.org

*atoms*:
    The list of atoms.
    
*timestep*:
    The time step.

*temperature*:
    The desired temperature, in Kelvin.

*taut*:
    Time constant for Berendsen temperature coupling.

*fixcm*:
    If True, the position and momentum of the center of mass is
    kept unperturbed.  Default: True.

*pressure*:
    The desired pressure, in bar (1 bar = 1e5 Pa).

*taup*:
    Time constant for Berendsen pressure coupling.

*compressibility*:
    The compressibility of the material, water 4.57E-5 bar-1, in bar-1


::

  # Room temperature simulation (300K, 0.1 fs time step, atmospheric pressure)
  dyn = NPTBerendsen(atoms, timestep=0.1*units.fs, temperature=300,
                   taut=0.1*1000*units.fs, pressure = 1.01325,
                   taup=1.0*1000*units.fs, compressibility=4.57e-5)

