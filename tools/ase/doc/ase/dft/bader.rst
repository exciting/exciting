==============
Bader Analysis
==============

Henkelman *et. al* have implemented a fast and robust algorithm for
calculating the electronic charges on individual atoms in molecules or
crystals, based on the Bader partitioning scheme [Bader]_. In that
method, the analysis is based purely on the electron density. The
partitioning of the density is determined according to its zero-flux
surfaces. Details of their implementation can be found in [Tang]_. The
program is freely available at
http://theory.cm.utexas.edu/henkelman/research/bader/ where you will
also find a description of the method.

The algorithm is very well suited for large solid state physical
systems, as well as large biomolecular systems. The computational time
depends only on the size of the 3D grid used to represent the electron
density, and scales linearly with the total number of grid points. As
the accuracy of the method depends on the grid spacing, it is
recommended to check for convergnce in this parameter (which should
usually by smaller than the default value).

The program takes cube input files. It does *not* support units, and
assumes atomic units for the density (Bohr^-3).

All ase dft calculators have a ``get_pseudo_density`` method, which
can be used to get the density. A simple python script for making a
cube file, ready for the Bader program, could be::

  >>> from ase import *
  >>> density = calc.get_pseudo_density() * Bohr**3
  >>> write('filename.cube', atoms, data=density)

Some calculators (e.g. gpaw) also have a method called
``get_all_electron_density``, in which case this is preferable to
``get_pseudo_density``.

Note that it is strongly recommended to use version 0.26b or higher of
the program, and the examples below refer to this version.

.. [Bader] R. F. W. Bader.  Atoms in Molecules: A Quantum Theory.
           Oxford University Press, New York, 1990.

.. [Tang]  W. Tang, E. Sanville, G. Henkelman.
           A grid-based Bader analysis algorithm without lattice bias.
           J. Phys.: Compute Mater. 21, 084204 (2009).


Example: The water molecule
---------------------------

The following example shows how to do Bader analysis for a water molecule.

First do a ground state calculation, and save the density as a cube file::

  from ase import *
  from gpaw import *

  atoms = molecule('H2O', cell=[7.5, 9, 9], calculator=GPAW(h=.17, xc='PBE'))
  atoms.center()
  atoms.get_potential_energy()

  rho = atoms.calc.get_all_electron_density(gridrefinement=4) * Bohr**3
  write('water_density.cube', atoms, data=rho)

Then analyse the density cube file by running (use `bader -h` for a
description of the possible options)::

  $ bader -p all_atom -p atom_index water_density.cube

This will produce a number of files. The `ACF.dat` file, contains a
summary of the Bader analysis::

  |     #         X           Y           Z        CHARGE     MIN DIST
  |  -----------------------------------------------------------------
  |     1      7.0865      8.5038      9.0672      9.1121      1.3250 
  |     2      7.0865      9.9461      7.9403      0.4440      0.2834 
  |     3      7.0865      7.0615      7.9403      0.4440      0.2834 
  |  -----------------------------------------------------------------
  |    NUMBER OF ELECTRONS:        9.99999

Revealing that 0.56 electrons have been transfered from each
Hydrogen atom to the Oxygen atom.

The `BvAtxxxx.dat` files, are cube files for each Bader volume,
describing the density within that volume. (I.e. it is just the
original cube file, truncated to zero outside the domain of the
specific Bader volume).

`AtIndex.dat` is a cube file with an integer value at each grid point,
describing which Bader volume it belongs to.

The plot below shows the dividing surfaces of the Hydrogen Bader
volumes. This was achieved by plotting a contour surface of
`AtIndex.dat` at an isovalue of 1.5.

.. image:: water_divide_surf.png
   :height: 220 pt

You can attach the output charges from the bader program to the atoms
for further processing::

  from ase import *
  from ase.io.bader import attach_charges

  # define the molecule as above
  atoms = molecule('H2O')
  atoms.set_cell([7.5, 9, 9])
  atoms.center()

  # the next two lines are equivalent (only one needed)
  attach_charges(atoms)
  attach_charges(atoms, 'ACF.dat')

  for atom in atoms:
      print 'Atom', atom.symbol, 'Bader charge', atom.charge
