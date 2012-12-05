.. module:: ase.atom

The Atom object
===============

ASE defines a python class called :class:`Atom` to setup and handle atoms
in electronic structure and molecular simulations. From a python
script, atoms can be created like this:

>>> from ase import Atom
>>> a1 = Atom('Si', (0, 0, 0))
>>> a2 = Atom('H', (1.3, 0, 0), mass=2)
>>> a3 = Atom(position=(0, 0, 0), Z=14)  # same is a1

.. autoclass:: Atom

The first argument to the constructor of an :class:`Atom` object is
the chemical symbol, and the second argument is the position in Å
units (see :mod:`units`).  The position can be any numerical sequence
of length three.  The properties of an atom can also be set using
keywords like it is done in the *a2* and *a3* examples above.

More examples:

>>> a = Atom('O', charge=-2)
>>> b = Atom(8, charge=-2)
>>> c = Atom('H', (1, 2, 3), magmom=1)
>>> print a.charge, a.position
-2 [ 0. 0. 0.]
>>> c.x = 0.0
>>> c.position
array([ 0.,  2.,  3.])
>>> b.symbol
'O'
>>> c.tag = 42
>>> c.number
1
>>> c.symbol = 'Li'
>>> c.number
3

If the atom object belongs to an Atoms object, then assigning
values to the atom attributes will change the corresponding
arrays of the atoms object:

>>> OH = Atoms('OH')
>>> OH[0].charge = -1
>>> OH.get_charges()
array([-1.,  0.])

Another example:

>>> for atom in bulk:
...     if atom.symbol = 'Ni':
...         atom.magmom = 0.7  # set initial magnetic moment
    

The different properties of an atom can be obtained and changed via
attributes (``position``, ``number``, ``tag``, ``momentum``, ``mass``,
``magmom``, ``charge``, ``x``, ``y``, ``z``):

>>> a1.position = [1, 0, 0]
>>> a1.position
array([ 1.,  0.,  0.])
>>> a1.z = 2.5
>>> a1.position
array([ 1. ,  0. ,  2.5])
>>> a2.magmom = 1.0

That last line will set the initial magnetic moment that some
calulators use (similar to the
:meth:`~ase.atoms.Atoms.set_initial_magnetic_moments` method).


.. note::

   The ``position`` and ``momentum`` attributes refer to mutable
   objects, so in some cases, you may want to use
   ``a1.position.copy()`` in order to avoid changing the position of
   ``a1`` by accident.



Getting an Atom from an Atoms object
------------------------------------

Indexing an :class:`Atoms` object returns an :class:`Atom` object
still remembering that it belongs to the collective :class:`Atoms`:
Modifying it will also change the atoms object:

>>> atoms = ase.data.molecules.molecule('CH4')
>>> atoms.get_positions()
array([[ 0.      ,  0.      ,  0.      ],
       [ 0.629118,  0.629118,  0.629118],
       [-0.629118, -0.629118,  0.629118],
       [ 0.629118, -0.629118, -0.629118],
       [-0.629118,  0.629118, -0.629118]])
>>> a = atoms[2]
>>> a
Atom('H', [-0.62911799999999996, -0.62911799999999996, 0.62911799999999996], index=2)
>>> a.x = 0
>>> atoms.get_positions()
array([[ 0.      ,  0.      ,  0.      ],
       [ 0.629118,  0.629118,  0.629118],
       [ 0.      , -0.629118,  0.629118],
       [ 0.629118, -0.629118, -0.629118],
       [-0.629118,  0.629118, -0.629118]])
                                                   

.. seealso::

   :epydoc:`atom.Atom`:
     All the details!

   :mod:`atoms`:
     More information about how to use collections of atoms.

   :mod:`calculators`:
     Information about how to calculate forces and energies of atoms.

