.. module:: calculators
   :synopsis: Energy, force and stress calculators.

===========
Calculators
===========

For ASE, a calculator is a black box that can take atomic numbers and
atomic positions from an :class:`~atoms.Atoms` object and calculate the
energy and forces and sometimes also stresses.

In order to calculate forces and energies, you need to attach a
calculator object to your atoms object:

>>> a = read('molecule.xyz')
>>> e = a.get_potential_energy()
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/jjmo/ase/ase/atoms.py", line 399, in get_potential_energy
    raise RuntimeError('Atoms object has no calculator.')
RuntimeError: Atoms object has no calculator.
>>> from ase.calculators import Abinit
>>> calc = Abinit(...)
>>> a.set_calculator(calc)
>>> e = a.get_potential_energy()
>>> print e
-42.0

Here, we used the :meth:`~ase.atoms.Atoms.set_calculator` method to attach
an instance of the :class:`~abinit.Abinit` class and then
we asked for the energy.

Alternatively, a calculator can be attached like this::

  atoms = Atoms(..., calculator=Siesta())


.. _supported calculators:

Supported calculators
=====================


================  ===========================================  ============
Code              Description                                  Type
================  ===========================================  ============
GPAW_             Grid-based real-space PAW code               :term:`DFT`,
                                                               :term:`HF`
Asap_             Highly efficient EMT code (written in C++)   :term:`EMT`
:mod:`jacapo`     ASE interface to Dacapo_,                    :term:`DFT`
                  a planewave ultra-soft pseudopotential code
Dacapo_           Old interface to Dacapo_. Requires           :term:`DFT`
                  Numeric python and ASE2.
:mod:`emt`        Effective Medium Theory calculator           :term:`EMT`
:mod:`abinit`     A planewave pseudopotential code             :term:`DFT`
:mod:`siesta`     LCAO pseudopotential code                    :term:`DFT`
:mod:`dftb`       DftbPlus_ DFT based tight binding            :term:`DFT`
:mod:`turbomole`  Fast atom orbital code Turbomole_,           :term:`DFT`,
                                                               :term:`HF`
:mod:`castep`     Planewave pseodopotential code               :term:`DFT`,
                                                               :term:`HF`
:mod:`vasp`       Planewave PAW code                           :term:`DFT`
:mod:`FHI-aims`   Numeric Atomic Orbital, full potential code  :term:`DFT`,
                                                               :term:`HF` 
:mod:`exciting`   Full Potential LAPW code                     :term:`DFT`,
                                                               :term:`LAPW`
:mod:`fleur`      Full Potential LAPW code                     :term:`DFT`,
                                                               :term:`LAPW`
:mod:`lammps`     Classical molecular dynamics code
:mod:`mmtk`       XXX Library for molecular simulations 
================  ===========================================  ============
  

.. _Asap: http://wiki.fysik.dtu.dk/asap
.. _Dacapo: http://wiki.fysik.dtu.dk/dacapo
.. _GPAW: http://wiki.fysik.dtu.dk/gpaw
.. _DftbPlus: http://www.dftb-plus.info/
.. _Turbomole: http://www.turbomole.com/


The calculators can be divided in three groups:

1) GPAW, Asap, Dacapo have their own native ASE interfaces.

2) Jacapo, ABINIT, SIESTA, DftbPlus, TURBOMOLE, VASP, FLEUR, 
   FHI-aims, LAMMPS and MMTK have Python wrappers in the ASE
   package, but the actual codes are not part of ASE.

3) EMT is a pure python implementation of the Effective Medium Theory
   potential and it is included in the ASE package.


Documentation for group 2 and 3 calculators
===========================================

.. toctree::

   emt
   jacapo
   abinit
   siesta
   dftb
   turbomole
   vasp
   FHI-aims
   lammps
   mmtk
   exciting
   fleur
   castep



Calculator interface
====================

All calculators must have the following interface:

.. autoclass:: ase.calculators.interface.Calculator
   :members:


Electronic structure calculators
================================

These calculators have wave functions, electron densities, eigenvalues
and many other quantities.  Therefore, it makes sense to have a set of
standard methods for accessing those quantities:

.. autoclass:: ase.calculators.interface.DFTCalculator
   :members:



Building new calculators
========================

Adding an ASE interface to your favorite force-calculator can be very
simple.  Take a look at the Python wrapper we have in the ASE code for
using the SIESTA_ code with ASE: :trac:`ase/calculators/siesta.py`.
Here, a :class:`~siesta.Siesta` class is defined.  An instance of this class
will simply write the fdf input-file, start the SIESTA Fortran
program, and finally read the energy, forces and stresses from the
text output-file.


.. _SIESTA: http://www.uam.es/departamentos/ciencias/fismateriac/siesta/



Building neighbor-lists
=======================

The :class:`EMT` potential and the GPAW_ DFT calculator both make
use of ASE's built-in neighbor-list class:

.. autoclass:: ase.calculators.neighborlist.NeighborList
   :members:


