.. _atommanip:

Manipulating atoms
------------------
**XXX rewrite this:  use fcc111() function and do some relaxation ...**

We will set up a one layer slab of Ni atoms with one Ag adatom.

Define the slab atoms:

>>> from ase import Atoms
>>> atoms = Atoms('Ni4', [(0, 0, 0),
...                       (0.45, 0, 0),
...                       (0, 0.5, 0),
...                       (0.5, 0.5, 0)])

Have a look at the individual atoms:

>>> atoms[0]
Atom('Ni', [0.0, 0.0, 0.0], atoms=..., index=0)
>>> atoms[1]
Atom('Ni', [0.45, 0.0, 0.0], atoms=..., index=1)
>>> atoms[2]
Atom('Ni', [0.0, 0.5, 0.0], atoms=..., index=2)
>>> atoms[3]
Atom('Ni', [0.5, 0.5, 0.0], atoms=..., index=3)

Let us assume we forgot how many atoms we set up:

>>> atoms[4]
Traceback (most recent call last):
File "<stdin>", line 1, in ?
IndexError: list index out of range

Wrong because we only have four atoms

>>> len(atoms)
4

Change the position of the 2nd atom in the list

>>> atoms[1].x = 0.5
>>> atoms.get_positions()
array([[ 0. ,  0. ,  0. ],
       [ 0.5,  0. ,  0. ],
       [ 0. ,  0.5,  0. ],
       [ 0.5,  0.5,  0. ]])

What is the unit cell so far?

>>> atoms.get_cell()
array([[ 1.,  0.,  0.],
       [ 0.,  1.,  0.],
       [ 0.,  0.,  1.]])

Now, setup a p(2x2) cell in a hexagonal surface.
Here, *a* is the fcc lattice constant, the cell is 10 layers high:

>>> from numpy import sqrt
>>> a = 3.55
>>> cell = [(2/sqrt(2.)*a, 0, 0),
...         (1/sqrt(2.)*a, sqrt(3./2.)*a, 0),
...         (0, 0, 10*sqrt(3.)/3.*a)]
>>> cell
[(5.0204581464244864, 0, 0),
(2.5102290732122432, 4.3478442934401409, 0),
(0, 0, 20.495934556231713)]
>>> atoms.set_cell(cell, scale_atoms=True)

The argument *scale_atoms=True* indicates that the atomic positions should be
scaled with the unit cell. The default is *scale_atoms=False* indicating that
the cartesian coordinates remain the same when the cell is changed.

>>> atoms.get_positions()
array([[ 0.        ,  0.        ,  0.        ],
       [ 2.51022907,  0.        ,  0.        ],
       [ 1.25511454,  2.17392215,  0.        ],
       [ 3.76534361,  2.17392215,  0.        ]])

Plot the whole system by bringing up the :mod:`gui`:

>>> from ase.visualize import view
>>> view(atoms)

.. image:: a1.png
   :scale: 35

Within the viewer (called :mod:`ag <gui>` or :mod:`ase.gui <gui>`) it
is possible to repeat the unit cell in all three directions (using the
:menuselection:`Repeat --> View` window).

.. image:: a2.png
   :scale: 35

We now add an adatom.  Since the supercell is now declared as the unit
cell for our atoms we can either add the atom using its cartesian
coordinates in Angstrom or rescale the unit cell and use scaled
coordinates. We try the latter:

>>> from numpy import identity
>>> from ase import Atom
>>> xyzcell = identity(3) # The 3x3 unit matrix
>>> atoms.set_cell(xyzcell, scale_atoms=True)  # Set the unit cell and rescale
>>> atoms.append(Atom('Ni', (1/6., 1/6., .1)))  
>>> atoms.set_cell(cell, scale_atoms=True)  # Set the unit cell and scale back

The structure now looks like this:

>>> view(atoms)

.. image:: a3.png
   :scale: 35
