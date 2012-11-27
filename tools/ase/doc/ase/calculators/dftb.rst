.. module:: dftb

=========
DftbPlus
=========

Introduction
============

DftbPlus_ is a density-functional based tight-binding code using 
atom centered orbitals. This 
interface makes it possible to use DftbPlus_ as a calculator in ASE.
You need to register at DftbPlus_ site to download the code.
Additionally you need Slater-coster files for the combination of 
atom types of your system. These can be obtained at dftb.org_.

.. _DftbPlus: http://www.dftb-plus.info/
.. _dftb.org: http://www.dftb.org/



Environment variables
=====================

Set environment variables in your configuration file (what is the directory
for the Slater-Koster files and what is the name of the executable):

- bash::

  $ DFTB_PREFIX=/my_disk/my_name/lib/Dftb+sk/mio-0-1/  (an example)
  $ DFTB_COMMAND=~/bin/DFTB+/dftb+_s081217.i686-linux  (an example)

- csh/tcsh::

  $ setenv DFTB_PREFIX /my_disk/my_name/lib/Dftb+sk/mio-0-1/  (an example)
  $ setenv DFTB_COMMAND ~/bin/DFTB+/dftb+_s081217.i686-linux   (an example)


DftbPlus Calculator
==================== 
This is a preliminary version of the DftbPlus calculator, so all the
DftbPlus features are unfortunately not there.

Use write_dftb=False to use your own 'dftb_in.hsd'-file with all the
flavours you need. In that case ase only updates the coordinates of
the system, otherwise file 'dftb_in.hsd' remains intact (see example 2
below). However, in case of own 'dftb_in.hsd' file you need first read
in atom position and atom type information of your system for ase (for
instance a xyz-file), see example 2 below.

The atom positions in file 'dftb_in.hsd' are updated during 
ASE geometry optimization.

For the spin polarised calculations this ASE-interface generates parameters 
with GGA-PBE-spin parameters. If you need LDA use your own 'dftb_in.hsd'-file. 

Information of periodicity is taken from ase (see example1 below).


For example::

    calc = Dftb(label='o2',
                write_dftb=True,
                do_spin_polarized=True,
                unpaired_electrons=2.0,
                fermi_temperature=100.0,
                scc=True)	


Parameters
==========
label: str
    Prefix to use for filenames (label.txt, ...).
    Default is 'dftb'.
write_dftb: boolean
    True: a minimal input file (name of which is always 'dftb_in.hsd')
    is written based on values given here.
    False: input file for dftb+ is not written. User must have
    generated file 'dftb_in.hsd' in the working directory.
    Use write_dftb=False to use your own 'dftb_in.hsd'-file.
charge: float
    Total charge of the system.
include_dispersion: boolean
    True: Default dispersion parameters are written in the 
    file 'dftb_in.hsd' (requires that also write_dftb_input_file==True)
    False: dispersion parameters are not written here.
do_spin_polarized: boolean
    True: Spin polarized calculation
    False: Spin unpolarized calculation
unpaired_electrons: float
    Number of spin unpaired electrons in the system.
    Relevant only if do_spin_polarized==True
fermi_temperature: float
    Fermi temperature for electrons.
scc: boolean
    True: Do charge self consistent dftb+
    False: No SCC, charges on atoms are not iterated
    
Example1: Geometry Optimization
===============================

.. literalinclude:: dftb_ex1_relax.py

Example2: Geometry Optimization using your own 'dftb_in.hsd' file
==================================================================

.. literalinclude:: dftb_ex2_relax.py

Input file for example 2 (h2o_1.xyz)

.. literalinclude:: h2o_1.xyz

You also need to have generated the file 'dftb_in.hsd', here is and example:
change the variable 'Prefix' below 

.. literalinclude:: dftb_in.hsd




