===================
Nudged elastic band
===================

.. module:: neb
   :synopsis: Nudged Elastic Band method.

.. default-role:: math

The Nudged Elastic Band method is a technique for finding transition paths
(and corresponding energy barriers) between given initial and final states.
The method involves constructing a "chain" of "replicas" or "images" of the
system and relaxing them in a certain way.

Relevant literature References:


1. H. Jonsson, G. Mills, and K. W. Jacobsen, in 'Classical and Quantum
   Dynamics in Condensed Phase Systems', edited by B. J. Berne,
   G. Cicotti, and D. F. Coker, World Scientific, 1998 [standard
   formulation]

2. 'Improved Tangent Estimate in the NEB method for Finding Minimum
   Energy Paths and Saddle Points', Graeme Henkelman and Hannes
   Jonsson, J. Chem. Phys. 113, 9978 (2000) [improved tangent
   estimates]

3. 'A Climbing-Image NEB Method for Finding Saddle Points and Minimum
   Energy Paths', Graeme Henkelman, Blas P. Uberuaga and Hannes
   Jonsson, J. Chem. Phys. 113, 9901 (2000)


The NEB class
=============

This module defines one class:

.. autoclass:: ase.neb.NEB

Example of use, between initial and final state which have been previously
saved in A.traj and B.traj::

  from ase import io
  from ase.neb import NEB
  from ase.optimize import MDMin
  # Read initial and final states:
  initial = io.read('A.traj')
  final = io.read('B.traj')
  # Make a band consisting of 5 images:
  images = [initial]
  images += [initial.copy() for i in range(3)]
  images += [final]
  neb = NEB(images)
  # Interpolate linearly the potisions of the three middle images:
  neb.interpolate()
  # Set calculators:
  for image in images[1:4]:
      image.set_calculator(MyCalculator(...))
  # Optimize:
  optimizer = MDMin(neb, trajectory='A2B.traj')
  optimizer.run(fmax=0.04)

Be sure to use the copy method (or similar) to create new instances
of atoms within the list of images fed to the NEB. Do *not* use something
like [initial for i in range(3)], as it will only create references to
the original atoms object.

Notice the use of the :meth:`~NEB.interpolate` method to get a good
initial guess for the path from A to B.

.. method:: NEB.interpolate()

   Interpolate path linearly from initial to final state.

Only the internal images (not the endpoints) need have
calculators attached.


.. seealso::

   :mod:`optimize`:
        Information about energy minimization (optimization).

   :mod:`calculators`:
        How to use calculators.

   :ref:`tutorials`:

        * :ref:`diffusion_tutorial`
        * :ref:`neb1`
        * :ref:`neb2`


.. note::

  If there are `M` images and each image has `N` atoms, then the NEB
  object behaves like one big Atoms object with `MN` atoms, so its
  :meth:`get_positions` method will return a `MN \times 3` array.


Trajectories
============

The code::

  from ase.optimize import QuasiNewton
  optimizer = QuasiNewton(neb, trajectory='A2B.traj')

will write all images to one file.  The Trajectory object knows about
NEB calculations, so it will write `M` images with `N` atoms at every
iteration and not one big configuration containing `MN` atoms.

The result of the latest iteration can now be analysed with this
command: :command:`ag A2B.traj@-5:`.

For the example above, you can write the images to individual
trajectory files like this::

  for i in range(1, 4):
      qn.attach(io.PickleTrajectory('A2B-%d.traj' % d, 'w', images[i])

The result of the latest iteration can be analysed like this:

.. highlight:: bash

::

  $ ag A.traj A2B-?.traj B.traj -n -1 

.. highlight:: python


Restarting
==========

Restart the calculation like this::

  images = io.read('A2B.traj@-5:')



Climbing image
==============

The "climbing image" variation involves designating a specific image to behave
differently to the rest of the chain: it feels no spring forces, and the
component of the potential force parallel to the chain is reversed, such that
it moves towards the saddle point. This depends on the adjacent images
providing a reasonably good approximation of the correct tangent at the
location of the climbing image; thus in general the climbing image is not
turned on until some iterations have been run without it (generally 20% to 50%
of the total number of iterations).

To use the climbing image NEB method, instantiate the NEB object like this::

  neb = NEB(images, climb=True)

.. note::

  Quasi-Newton methods, such as BFGS, are not well suited for climbing image
  NEB calculations. FIRE have been known to give good results, although
  convergence is slow.


Parallelization over images
===========================

Some calculators can parallelize over the images of a NEB calculation.
The script will have to be run with an MPI-enabled Python interpreter
like GPAW_'s gpaw-python_.  All images exist on all processors, but
only some of them have a calculator attached::

  from ase.parallel import rank, size
  from ase.calculators.emt import EMT
  # Number of internal images:
  n = len(images) - 2
  j = rank * n // size
  for i, image in enumerate(images[1:-1]):
      if i == j:
          image.set_calculator(EMT())

Create the NEB object with ``NEB(images, parallel=True)`` and let the
master processes write the images::

  if rank % (size // n) == 0:
      traj = io.PickleTrajectory('neb%d.traj' % j, 'w', images[1 + j],
                                 master=True)
      optimizer.attach(traj)

For a complete example using GPAW_, see here_.

.. _GPAW: http://wiki.fysik.dtu.dk/gpaw
.. _gpaw-python: https://wiki.fysik.dtu.dk/gpaw/documentation/manual.html#parallel-calculations
.. _here: https://wiki.fysik.dtu.dk/gpaw/tutorials/neb/neb.html

.. default-role::
