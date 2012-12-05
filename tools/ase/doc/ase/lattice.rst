.. default-role:: math

.. _general-crystal-section:

General crystal structures and surfaces
=======================================

.. module:: lattice

Modules for creating crystal structures are found in the module
:mod:`lattice`.  Most Bravais lattices are implemented, as
are a few important lattices with a basis.  The modules can create
lattices with any orientation (see below).  These modules can be used
to create surfaces with any crystal structure and any orientation by
later adding a vacuum layer with :func:`lattice.surface.add_vacuum`.

Example
-------

To set up a slab of FCC copper with the [1,-1,0] direction along the
x-axis, [1,1,-2] along the y-axis and [1,1,1] along the z-axis, use::

  from ase.lattice.cubic import FaceCenteredCubic
  atoms = FaceCenteredCubic(directions=[[1,-1,0], [1,1,-2], [1,1,1]],
                            size=(2,2,3), symbol='Cu', pbc=(1,1,0))

The minimal unit cell is repeated 2*2*3 times.  The lattice constant
is taken from the database of lattice constants in :mod:`data` module.
There are periodic boundary conditions along the *x* and *y* axis, but
free boundary conditions along the *z* axis. Since the three directions
are perpendicular, a (111) surface is created.

To set up a slab of BCC copper with [100] along the first axis, [010]
along the second axis, and [111] along the third axis use::

  from ase.lattice.cubic import BodyCenteredCubic
  atoms = BodyCenteredCubic(directions=[[1,0,0], [0,1,0], [1,1,1]],
                            size=(2,2,3), symbol='Cu', pbc=(1,1,0),
			    latticeconstant=4.0)

Since BCC is not the natural crystal structure for Cu, a lattice
constant has to be specified.  Note that since the repeat directions
of the unit cell are not orthogonal, the Miller indices of the
surfaces will *not* be the same as the Miller indices of the axes.
The indices of the surfaces in this example will be (1,0,-1), (0,1,-1)
and (0,0,1).


Available crystal lattices
--------------------------

The following modules are currently available (the * mark lattices
with a basis):

* ``lattice.cubic``

  - ``SimpleCubic`` 
  - ``FaceCenteredCubic``
  - ``BodyCenteredCubic``
  - ``Diamond`` (*)

* ``lattice.tetragonal``

  - ``SimpleTetragonal``
  - ``CenteredTetragonal``

* ``lattice.orthorhombic``

  - ``SimpleOrthorhombic``
  - ``BaseCenteredOrthorhombic``
  - ``FaceCenteredOrthorhombic``
  - ``BodyCenteredOrthorhombic``

* ``lattice.monoclinic``

  - ``SimpleMonoclinic``
  - ``BaseCenteredMonoclinic``

* ``lattice.triclinic``

  - ``Triclinic``

* ``lattice.hexagonal``

  - ``Hexagonal``
  - ``HexagonalClosedPacked`` (*)
  - ``Graphite`` (*)

* The rhombohedral (or trigonal) lattices are not implemented.  They
  will be implemented when the need arises (and if somebody can tell
  me_ the precise definition of the 4-number Miller indices - I only
  know that they are "almost the same as in hexagonal lattices").

* ``lattice.compounds``

  Lattices with more than one element.  These are mainly intended as
  examples allowing you to define new such lattices.  Currenly, the
  following are defined

  - ``B1`` = ``NaCl`` = ``Rocksalt``
  - ``B2`` = ``CsCl``
  - ``B3`` = ``ZnS`` = ``Zincblende``
  - ``L1_2`` = ``AuCu3``
  - ``L1_0`` = ``AuCu``

.. _me: http://www.fysik.dtu.dk/~schiotz

Usage
-----

The lattice objects are called with a number of arguments specifying
e.g. the size and orientation of the lattice.  All arguments should be
given as named arguments.  At a minimum the ``symbol`` argument must
be specified.


``symbol``
  The element, specified by the atomic number (an integer) or by the
  atomic symbol (i.e. 'Au').  For compounds, a tuple or list of
  elements should be given.  This argument is mandatory.

``directions`` and/or ``miller``: 
  Specifies the orientation of the
  lattice as the Miller indices of the three basis vectors of the
  supercell (``directions=...``) and/or as the Miller indices of the
  three surfaces (``miller=...``).  Normally, one will specify either
  three directions or three surfaces, but any combination that is both
  complete and consistent is allowed, e.g. two directions and two
  surface miller indices (this example is slightly redundant, and
  consistency will be checked).  If only some directions/miller
  indices are specified, the remaining should be given as ``None``.
  If you intend to generate a specific surface, and prefer to specify
  the miller indices of the unit cell basis (``directions=...``), it
  is a good idea to give the desired Miller index of the surface as
  well to allow the module to test for consistency.  Example:

  >>> atoms = BodyCenteredCubic(directions=[[1,-1,0],[1,1,-1],[0,0,1]],
  ...                           miller=[None, None, [1,1,2]], ...)

  If neither ``directions`` nor ``miller`` are specified, the default
  is ``directions=[[1,0,0], [0,1,0], [0,0,1]]``.

``size``:
  A tuple of three numbers, defining how many times the fundamental
  repeat unit is repeated. Default: (1,1,1).  Be aware that if
  high-index directions are specified, the fundamental repeat unit may
  be large.

``latticeconstant``:
  The lattice constant.  If no lattice constant is
  specified, one is extracted from ASE.ChemicalElements provided that
  the element actually has the crystal structure you are creating.
  Depending on the crystal structure, there will be more than one
  lattice constant, and they are specified by giving a dictionary or a
  tuple (a scalar for cubic lattices).  Distances are given in
  Angstrom, angles in degrees. 

  =============  ===================  ========================================
  Structure      Lattice constants    Dictionary-keys
  =============  ===================  ========================================
  Cubic          a                    'a'
  Tetragonal     (a, c)               'a', 'c' or 'c/a'
  Orthorhombic   (a, b, c)            'a', 'b' or 'b/a', 'c' or 'c/a'
  Triclinic      (a, b, c, `\alpha`,  'a', 'b' or 'b/a', 'c' or
                 `\beta`, `\gamma`)   'c/a', 'alpha', 'beta', 'gamma'
  Monoclinic     (a, b, c, alpha)     'a', 'b' or 'b/a', 'c' or 'c/a', 'alpha'
  Hexagonal      (a, c)               'a', 'c' or 'c/a'
  =============  ===================  ========================================
  
  Example:

  >>> atoms = Monoclinic( ... , latticeconstant={'a': 3.06, 
  ...     'b/a': 0.95, 'c/a': 1.07, 'alpha': 74})


``debug``:
  Controls the amount of information printed.  0: no info is printed.
  1 (the default): The indices of surfaces and unit cell vectors are
  printed.  2: Debugging info is printed.


Defining new lattices
---------------------

Often, there is a need for new lattices - either because an element
crystallizes in a lattice that is not a simple Bravais lattice, or
because you need to work with a compound or an ordered alloy.

All the lattice generating objects are instances of a class, you
generate new lattices by deriving a new class and instantiating it.
This is best explained by an example.  The diamond lattice is two
interlacing FCC lattices, so it can be seen as a face-centered cubic
lattice with a two-atom basis.  The Diamond object could be defined like
this::

  from ase.lattice.cubic import FaceCenteredCubicFactory
  class DiamondFactory(FaceCenteredCubicFactory):
      """A factory for creating diamond lattices."""
      xtal_name = 'diamond'
      bravais_basis = [[0, 0, 0], [0.25, 0.25, 0.25]]
    
  Diamond = DiamondFactory()



Lattices with more than one element
```````````````````````````````````

Lattices with more than one element is made in the same way.  A new
attribute, ``element_basis``, is added, giving which atoms in the
basis are which element.  If there are four atoms in the basis, and
element_basis is (0,0,1,0), then the first, second and fourth atoms
are one element, and the third is the other element.  As an example,
the AuCu3 structure (also known as `\mathrm{L}1_2`) is defined as::

  # The L1_2 structure is "based on FCC", but is really simple cubic
  # with a basis.
  class AuCu3Factory(SimpleCubicFactory):
      "A factory for creating AuCu3 (L1_2) lattices."
      bravais_basis = [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
      element_basis = (0, 1, 1, 1)

  AuCu3 = L1_2 = AuCu3Factory()

Sometimes, more than one crystal structure can be used to define the
crystal structure, for example the Rocksalt structure is two
interpenetrating FCC lattices, one with one kind of atoms and one with
another.  It would be tempting to define it as

::

  class NaClFactory(FaceCenteredCubicFactory):
      "A factory for creating NaCl (B1, Rocksalt) lattices."

      bravais_basis = [[0, 0, 0], [0.5, 0.5, 0.5]]
      element_basis = (0, 1)


  B1 = NaCl = Rocksalt = NaClFactory()

but if this is used to define a finite system, one surface would be
covered with one type of atoms, and the opposite surface with the
other.  To maintain the stochiometry of the surfaces, it is better to
use the simple cubic lattice with a larger basis::

  # To prevent a layer of element one on one side, and a layer of
  # element two on the other side, NaCl is based on SimpleCubic instead
  # of on FaceCenteredCubic
  class NaClFactory(SimpleCubicFactory):
      "A factory for creating NaCl (B1, Rocksalt) lattices."

      bravais_basis = [[0, 0, 0], [0, 0, 0.5], [0, 0.5, 0], [0, 0.5, 0.5],
		       [0.5, 0, 0], [0.5, 0, 0.5], [0.5, 0.5, 0],
		       [0.5, 0.5, 0.5]]
      element_basis = (0, 1, 1, 0, 1, 0, 0, 1)


  B1 = NaCl = Rocksalt = NaClFactory()

More examples can be found in the file :trac:`ase/lattice/compounds.py`.

.. default-role::
