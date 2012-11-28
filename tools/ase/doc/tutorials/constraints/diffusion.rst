.. _constraints_diffusion_tutorial:

======================================================
Diffusion of gold atom on Al(100) surface (constraint)
======================================================

In this tutorial, we will calculate the energy barrier that was found
using the :mod:`NEB <neb>` method in the :ref:`diffusion_tutorial`
tutorial.  Here, we use a siple :class:`~ase.constraints.FixedPlane`
constraint that forces the Au atom to relax in the *yz*-plane only:

.. literalinclude:: diffusion4.py

The result can be analysed with the command :command:`ag mep?.traj -n
-1` (choose :menuselection:`Tools --> NEB`).  The barrier is found to
be 0.35 eV - exactly as in the :ref:`NEB <diffusion_tutorial>`
tutorial.

Here is a side-view of the path (unit cell repeated twice):

.. image:: diffusion-path.png


.. seealso::

   * :mod:`neb`
   * :mod:`constraints`
   * :ref:`diffusion_tutorial`
   * :func:`~lattice.surface.fcc100`
