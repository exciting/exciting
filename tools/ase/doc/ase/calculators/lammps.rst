.. module:: lammps

======
LAMMPS
======

Introduction
============

LAMMPS_ (Large-scale Atomic/Molecular Massively Parallel Simulator) is a classical molecular dynamics code.

    "LAMMPS has potentials for soft materials (biomolecules, polymers) and solid-state materials (metals, semiconductors) and coarse-grained or mesoscopic systems. It can be used to model atoms or, more generically, as a parallel particle simulator at the atomic, meso, or continuum scale."


.. _LAMMPS: http://lammps.sandia.gov



Environment variables
=====================

The environment variable :envvar:`LAMMPS_COMMAND` should contain
the path to the lammps binary, or more generally, a command line 
possibly also including an MPI-launcher command.
For example (in a Bourne-shell compatible environment):

.. highlight:: bash
 
::

  $ export LAMMPS_COMMAND=/path/to/lmp_binary

.. highlight:: python

or possibly something similar to

.. highlight:: bash
 
::

  $ export LAMMPS_COMMAND="/path/to/mpirun --np 4 lmp_binary"

.. highlight:: python



LAMMPS Calculator
================= 

The LAMMPS calculator first appeared in ASE version 3.5.0. 
At the time of the release of ASE 3.5.0, the LAMMPS calculator 
is still in a fairly early stage of development
(if you are missing some feature, consider checking out 
the source code development tree or some more recent version of ASE).

.. class:: LAMMPS(..., parameters={}, files=[], ...)

Below follows a list with a selection of parameters

==============  =========  ==============  =============================
keyword         type       default value   description
==============  =========  ==============  =============================
``files``       ``list``   ``[]``          List of files needed by 
                                           LAMMPS. Typically a list of
                                           potential files.
``parameters``  ``dict``   ``{}``          Dictionary with key-value
                                           pairs corresponding to 
                                           commands and arguments.
                                           Command-argument pairs 
					   provided here will
                                           be used for overriding the
                                           the calculator defaults.
==============  =========  ==============  =============================


Example
=======

A simple example.

::

  from ase import Atoms, Atom
  from ase.calculators.lammps import LAMMPS
  
  a = [6.5, 6.5, 7.7]
  d = 2.3608
  NaCl = Atoms([Atom('Na', [0, 0, 0]),
                Atom('Cl', [0, 0, d])],
               cell=a, pbc=True)
  
  calc = LAMMPS()
  NaCl.set_calculator(calc)
  
  print NaCl.get_stress()




