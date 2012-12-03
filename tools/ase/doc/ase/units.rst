.. module:: units

=====
Units
=====

Physical units are defined in the :trac:`ase/units.py` module.  Electron volts
(``eV``) and angstroms (``Ang``) are defined as 1.0.
Other units are
``nm``, ``Bohr``, ``Hartree`` or ``Ha``, ``kJ``, ``kcal``, ``mol``,
``Rydberg`` or ``Ry``, ``second``, ``fs`` and ``kB``.

.. note::

    All constants are taken from the 1986 CODATA_.

.. _CODATA: http://physics.nist.gov/cuu/Constants/archive1986.html

Examples:

>>> from ase import *
>>> 2 * Bohr
1.0583545150138329
>>> 25 * Rydberg
340.14244569396635
>>> 100 * kJ/mol
1.0364272141304978
>>> 300 * kB
0.025852157076770025
>>> 0.1 * fs
0.009822693531550318
>>> print '1 Hartree = '+str(Hartree*mol/kcal)+' kcal/mol'

=======================
The ``ase.data`` module
=======================

This module defines the following variables: ``atomic_masses``,
``atomic_names``, ``chemical_symbols``, ``covalent_radii``,
``cpk_colors`` and ``reference_states``.  All of these are lists that
should be indexed with an atomic number:

>>> atomic_names[92]
'Uranium'
>>> atomic_masses[2]
4.0026000000000002

If you don't know the atomic number of some element, then you can look
it up in the ``atomic_numbers`` dictionary:

>>> atomic_numbers['Cu']
29
>>> covalent_radii[29]
1.1699999999999999
