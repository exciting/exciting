.. module:: visualize

Visualization
=============

.. function:: view(atoms, data=None, viewer=None, repeat=None)

This provides an interface to various visualization tools, such as
:mod:`ase.gui <gui>`, :mod:`ase.visualize.vtk <vtk>`, RasMol_, VMD_, VTK_, gOpenMol_, or
Avogadro_. The default viewer is the ase.gui, described in the
:mod:`gui` module. The simplest invocation is::

  >>> from ase import view
  >>> view(atoms)

where ``atoms`` is any :class:`Atoms` object.  Alternative viewers can
be used by specifying the optional keyword ``viewer=...`` - use one of
'ase.gui', 'gopenmol', 'vmd', or 'rasmol'.  The VMD and Avogadro
viewers can take an optional ``data`` argument to show 3D data, such
as charge density::

  >>> view(atoms, viewer='VMD', data=array)

If you do not wish to open an interactive gui, but rather visualize
your structure by dumping directly to a graphics file; you can use the
``write`` command of the :mod:`io` module, which can write 'eps',
'png', and 'pov' files directly, like this::

  >>> write('image.png', atoms)

.. _RasMol: http://openrasmol.org/
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _VTK: http://www.vtk.org/VTK/project/about.html
.. _gOpenMol: http://www.csc.fi/gopenmol/
.. _Avogadro: http://avogadro.openmolecules.net/


VTK
------------

.. module:: visualize.vtk

The Visualization Toolkit (VTK_) is a powerful platform-independent graphics
engine, which comes as an open source graphics toolkit licensed under the 
`BSD license`_. It is available for a wide range of programming languages, 
including easily scriptable interfaces in Python and Tcl_.

.. _BSD license: http://en.wikipedia.org/wiki/BSD_licenses
.. _Tcl: http://www.tcl.tk/about

In the scientific community, VTK is used by thousands of researchers and
developers for 3D computer graphics, image processing, and visualization.
VTK includes a suite of 3D interaction widgets within the development
framework for information visualization, integrating GUI toolkits such as
Qt_ and Tk_ into a highly flexible design platform.

For visualization purposes within ASE, two different VTK-approaches are
supported, namely:

:Scripted on-the-fly rendering:
	ASE includes VTK-scripting for easy data visualization using the
	:mod:`vtk` module. Development is in progress, so you might want to
	check out the latest development release from SVN 
	(see :ref:`latest_development_release`).

:Interactive rendering:
	MayaVi_ is an easy-to-use GUI for VTK. With Enthought's traits-based
	VTK-wrapper (TVTK_), constructing VTK pipelines has been simplified greatly
	by introducing three basic concepts: data sources, filters and visualization
	modules. MayaVi also supports the VTK file formats, including the flexible
	VTK XML, which in ASE can be used to export atomic positions, forces and 
	volume data using the ``write`` command in the :mod:`io` module.

.. XXX -	`VTK Designer`_ is a visual editor for creating and editing VTK pipelines.

.. _Qt: http://www.qtsoftware.com/products
.. _Tk: http://www.tcl.tk/about
.. _VTK Designer: http://www.vcreatelogic.com/oss/vtkdesigner
.. _MayaVi: http://code.enthought.com/projects/mayavi
.. _TVTK: https://svn.enthought.com/enthought/wiki/TVTK

A key feature of VTK is the inherent ability to use MPI_ for parallel rending,
which is provided with built-in parallel composite rendering objects to handle
domain decomposition and subsequent recombination of the raster information.
This is particularly useful for non-interactive ray tracing, batch isosurface
generation and in-situ visualization of simulation data in cluster computing.

.. seealso::
	ParaView_ is a VTK-based open-source, multi-platform data analysis and
	visualization application for extremely large data-sets using distributed 
	memory computing resources and parallel rendering through MPI_.

.. _MPI: http://www.mpi-forum.org
.. _ParaView: http://www.paraview.org



PrimiPlotter
------------

The PrimiPlotter is intended to do on-the-fly plotting of the
positions of the atoms during long molecular dynamics simulations.
The module :mod:`ase.visualize.primiplotter` contains the PrimiPlotter
and the various output modules, see below.


.. autoclass:: ase.visualize.primiplotter.PrimiPlotter
   :inherited-members:



FieldPlotter
------------

The FieldPlotter is intended to plot fields defined on the atoms in
large-scale simulations.  The fields could be e.g. pressure, stress or
temperature (kinetic energy), i.e. any quantity that in a given
simulation is best defined on a per-atom basis, but is best
interpreted as a continuum field.

The current version of FieldPlotter only works if the number of atoms
is at least 5-10 times larger than the number of pixels in the plot.

.. autoclass:: ase.visualize.fieldplotter.FieldPlotter
   :inherited-members:




