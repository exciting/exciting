.. module:: io
   :synopsis: File input-output module


File input and output
=====================

The :mod:`io` module has two basic functions: :func:`read` and :func:`write`.
The two methods are described here:

.. autofunction:: ase.io.read
.. autofunction:: ase.io.write

The :func:`read` function is only designed to retrive the atomic configuration
from a file, but for the CUBE format you can import the function:

.. function:: read_cube_data


which will return a ``(data, atoms)`` tuple::
   
  from ase.io.cube import read_cube_data
  data, atoms = read_cube_data('abc.cube')



Examples
========

::

    from ase.lattice.surface import *
    adsorbate = Atoms('CO')
    adsorbate[1].z = 1.1
    a = 3.61
    slab = fcc111('Cu', (2, 2, 3), a=a, vacuum=7.0)
    add_adsorbate(slab, adsorbate, 1.8, 'ontop')
 
Write PNG image::

    write('slab.png', slab * (3, 3, 1), rotation='10z,-80x')

.. image:: io1.png

Write POVRAY file::

    write('slab.pov', slab * (3, 3, 1), rotation='10z,-80x')

This will write both a ``slab.pov`` and a ``slab.ini`` file.  Convert
to PNG with the command ``povray slab.ini`` or use the
``run_povray=True`` option:

.. image:: io2.png

Here is an example using ``bbox``::

    d = a / 2**0.5
    write('slab.pov', slab * (2, 2, 1),
          bbox=(d, 0, 3 * d, d * 3**0.5))

.. image:: io3.png

Note that the XYZ-format does not contain information about the unic cell:

>>> write('slab.xyz', slab)
>>> a = read('slab.xyz')
>>> a.get_cell()
array([[ 1.,  0.,  0.],
       [ 0.,  1.,  0.],
       [ 0.,  0.,  1.]])
>>> a.get_pbc()
array([False, False, False], dtype=bool)

Use ASE's native format for writing all information:

>>> write('slab.traj', slab)
>>> b = read('slab.traj')
>>> b.get_cell()
array([[  5.10531096e+00,  -4.11836034e-16,   1.99569088e-16],
       [  2.55265548e+00,   4.42132899e+00,   7.11236625e-17],
       [  8.11559027e+00,   4.68553823e+00,   1.32527034e+01]])
>>> b.get_pbc()
array([ True,  True,  True], dtype=bool)

A script showing all of the povray parameters, and generating the image below, 
can be found here: :trac:`doc/ase/save_pov.py`

.. image:: NaCl_C6H6.png

An other example showing how to change colors and textures in pov can
be found here: :trac:`doc/tutorials/saving_graphics.py`
