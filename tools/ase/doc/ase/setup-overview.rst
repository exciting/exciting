==========================
Creating atomic structures
==========================

.. toctree::
   :maxdepth: 2

   structure
   surface
   lattice
   molecules

ASE contains a number of modules for setting up atomic structures,
mainly molecules, bulk crystals and surfaces.  Some of these modules
have overlapping functionality, but strike a different balance between
flexibility and ease-of-use.


**Common bulk crystals**

The :func:`ase.structure.bulk` function can be used to create the most
common bulk crystal structures.  The function creates a single unit cell
oriented such that the number of atoms in the cell is minimal.

Read more: :ref:`bulk-crystal-section`.


**Common surfaces**

The :mod:`lattice.surface` module contains a number of 
functions for creating the most common surfaces in a minimal unit
cell, and for adding adsorbates to these surfaces.

Read more: :ref:`lattice-surface-section`.


**Nanotubes and nanoribbons**

The functions :func:`ase.structure.nanotube` and
:func:`ase.structure.graphene_nanoribbon` can be used to create Carbon
nanotubes and graphene sheets or nanoribbons.  Per default, they
create Carbon nanotubes and sheets, but other elements can be used. 

Read more:  :ref:`nanotubes-section` and :ref:`nanoribbons-section`.

**Generally oriented bulk crystals and surfaces**

The :mod:`lattice` module contains functions for creating most common
crystal structures with arbitrary orientation.  The user can specify
the desired Miller index along the three axes of the simulation, and
the smallest periodic structure fulfilling this specification is
created.  Thirteen of the 14 Bravais lattices are supported by the
module, as are a few lattices with a basis, and lattices for some of
the most common compounds/alloys.  The modules makes it possible to
define further lattices based on the supported Bravais lattices.

Both bulk crystals and surfaces can be created.

Read more: :ref:`general-crystal-section`.

**Molecules**

Some common molecules are available in the :mod:`~data.molecules` module.

Read more: :ref:`molecules-section`.
