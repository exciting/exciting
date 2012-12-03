.. _faq:

==========================
Frequently Asked Questions
==========================


ASE-GUI
=======

See also the :mod:`documantation for ag <gui>`.

How do I export images from a trajectory to png or pov files?
-------------------------------------------------------------

With ag, you can choose :menuselection:`File --> Save`, but this is
not fun if you need to do it for many images.  Here is how to do it on
the command line for a number of images::

  ag images.traj@0 -o image0.pov 
  ag images.traj@1 -o image1.pov 
  ag images.traj@2 -o image2.pov 

If you have many images, it will be easier to do it using the Python
interpreter:

>>> from ase import *
>>> for n, image in enumerate(read('images.traj@:3')):
...     write('image%d.pov' % n, image, run_povray=True, pause=False,
...           rotation='-90x,10z')

Here, we also:

* run povray to generate png files

* disable pausing between the images

* set a rotation (choose :menuselection:`View --> Rotate ...` in ag to select
  the best rotation angles)

Try:

>>> help(write)

to see all possibilities or read more :func:`here <ase.io.write>`.



General
=======

Citation: how should I cite ASE?
--------------------------------

If you find ASE useful in your research please cite:

   | S. R. Bahn and K. W. Jacobsen
   | `An object-oriented scripting interface to a legacy electronic structure code`__
   | Comput. Sci. Eng., Vol. **4**, 56-66, 2002

   __ http://dx.doi.org/10.1109/5992.998641

BibTex (:svn:`doc/ASE.bib`):

.. literalinclude:: ASE.bib


Download
========

Trying to checkout the code via SVN resulted::

 [~]$ svn checkout "https://svn.fysik.dtu.dk/projects/ase/trunk"
 svn: Unrecognized URL scheme 'https://svn.fysik.dtu.dk/projects/ase/trunk'

This error is diplayed in case the library 'libsvn_ra_dav' is missing on your system. The library is used by SVN, but is not installed by default. 
