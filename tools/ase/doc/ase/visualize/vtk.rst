.. module:: vtk
   :synopsis: Use of the Visualization Toolkit.

.. index:: vtk

===========
ASE-VTK
===========

For ASE, the :mod:`~visualize.vtk` interface consists of Python modules for 
automatic visualization of positions, bonds, forces and volume data (e.g. wave
functions) from an :class:`~atoms.Atoms` object, provided such data is made
available by the calculator.

.. note::

	The Python modules in ASE are intended to wrap lower-level functionality
	of the VTK object models in small and easy-to-comprehend classes. To be able
	to distinguish between build-in VTK objects and their wrappers, and because
	VTK uses the CamelCase naming convention whereas ASE uses lower-case cf. 
	our :ref:`coding conventions <python_codingstandard>`, all variables
	referring to VTK built-in types are prefixed by ``vtk_``. However, both VTK
	and wrapper classes are named according to the standard ``vtkFooBar``.


Representing atoms
==================
.. autoclass:: ase.visualize.vtk.atoms.vtkAtoms
   :members:
   :show-inheritance:

Atom-centered data
------------------

The superclass :class:`~ase.visualize.vtk.grid.vtkAtomicPositions` implements
the basic concepts for representing atomic-centered data in VTK.

.. autoclass:: ase.visualize.vtk.grid.vtkAtomicPositions
   :members:
..   :show-inheritance:

Predefined shapes
------------------

The class :class:`~ase.visualize.vtk.module.vtkGlyphModule` implements
the lower-level objects for representing predefined shapes (glyphs) in VTK.

.. autoclass:: ase.visualize.vtk.module.vtkGlyphModule
   :members:
   :inherited-members:
..   :show-inheritance:

