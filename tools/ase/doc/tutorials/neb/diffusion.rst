.. _diffusion_tutorial:

===============================================
Diffusion of gold atom on Al(100) surface (NEB)
===============================================

First, set up the initial and final states:

|initial| |final|

.. literalinclude:: diffusion1.py

.. note::  Notice how the tags are used to select the constrained atoms

Now, do the NEB calculation:

.. literalinclude:: diffusion2.py

|ts| |barrier|

.. note::

   For this reaction, the reaction coordinate is very simple: The
   *x*-coordinate of the Au atom.  In such cases, the NEB method is
   overkill, and a simple constraint method should be used like in this
   tutorial: :ref:`constraints_diffusion_tutorial`.

.. seealso::

   * :mod:`neb`
   * :mod:`constraints`
   * :ref:`constraints_diffusion_tutorial`
   * :func:`~lattice.surface.fcc100`
   


.. |initial| image:: diffusion-I.png
.. |final| image:: diffusion-F.png
.. |ts| image:: diffusion-T.png
.. |barrier| image:: diffusion-barrier.png


Parallelizing over images
=========================

Instead of having one process do the calculations for all three
internal images in turn, it will be faster to have three processes do
one image each.  This can be done like this:

.. literalinclude:: diffusion3.py


