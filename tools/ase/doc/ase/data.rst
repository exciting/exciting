.. module:: data

===============
The data module
===============


Atomic data
===========

This module defines the following variables:

.. data:: atomic_masses
.. data:: atomic_names
.. data:: chemical_symbols
.. data:: covalent_radii
.. data:: cpk_colors
.. data:: reference_states
.. data:: vdw_radii

All of these are lists that should be indexed with an atomic number:

>>> from ase.data import *
>>> atomic_names[92]
'Uranium'
>>> atomic_masses[2]
4.0026000000000002


.. data:: atomic_numbers

If you don't know the atomic number of some element, then you can look
it up in the :data:`atomic_numbers` dictionary:

>>> atomic_numbers['Cu']
29
>>> covalent_radii[29]
1.1699999999999999

The covalent radii are taken from [Cordeo08]_. 
The source of the van der Waals radii is given in vdw.py_.

.. [Cordeo08] *Covalent radii revisited*,
    Beatriz Cordero, Verónica Gómez, Ana E. Platero-Prats, Marc Revés,
    Jorge Echeverría, Eduard Cremades, Flavia Barragán and Santiago Alvarez,
    Dalton Trans., 2008, 2832-2838 DOI:10.1039/B801115J 

.. _vdw.py: https://trac.fysik.dtu.dk/projects/ase/browser/trunk/ase/data/vdw.py

.. _molecular data:

Molecular data
==============

The G1, G2, and G3-databases are available in the :mod:`molecules` module.

Example:

>>> from ase.data.molecules import molecule
>>> atoms = molecule('H2O')

All molecular members of each database is conveniently contained in a list
of strings (g1, g2, g3), and one can look up the
experimental atomization energy for each molecule.
This is extrapolated from experimental heats of formation at room temperature,
using calculated zero-point energies and thermal corrections.

Example:

>>> from ase.data.molecules import molecule, g2, get_atomization_energy
>>> g2[11]
'H2O'
>>> atoms = molecule(g2[11]) 
>>> get_atomization_energy(g2[11])
232.57990000000001
>>> from ase.units import kcal,mol
>>> get_atomization_energy(g2[11])*kcal/mol
10.08562144637833

where the last line converts the experimental atomization energy of H2O
from units of kcal/mol to eV.


S22, s26, and s22x5 data
========================

The s22, s26, and s22x5 databases are available in the :mod:`s22` module.

Each weakly bonded complex is identified as an entry in a list of strings
(s22, s26, s22x5), and is fully created by a 'create'-function:

>>> from ase.data.s22 import s22, create_s22_system
>>> sys = s22[0]
>>> sys
'Ammonia_dimer'
>>> atoms = create_s22_system(sys)
>>> atoms.get_chemical_symbols()
['N', 'H', 'H', 'H', 'N', 'H', 'H', 'H']

The coupled-cluster interaction energies for the s22 and s26 systems
are retrieved like this:

>>> from ase.data.s22 import s22, get_interaction_energy_s22
>>> get_interaction_energy_s22(s22[0])
-0.1375

in units of eV. For s22 these are not the original energies,
but from more recent work where the same (large) basis set
was used for all complexes, yielding more accurate
coupled-cluster interaction energies.

The s22x5 database expands on the original s22 data by introducing
non-equilibrium geometries for each complex
(0.9, 1.0, 1.2, 1.5, and 2.0 times original intermolecular distance).
However, these calculations were done in accordance with the methods
used in the original s22 work, and so is expected to inherit the
same problems with mixed basis set sizes.
Assuming the interaction energy error due to this is the same in all
5 geometries for each complex, the default s22x5 interaction energies
are therefore corrected with the energy difference between
original and newer energies at the original separation.

Example:

>>> from ase.data.s22 import *
>>> sys1 = s22[0]
>>> sys1
'Ammonia_dimer'
>>> atoms1 = create_s22_system(sys1)
>>> sys2 = s22x5[0]
>>> sys2
'Ammonia_dimer_0.9'
>>> atoms2 = create_s22_system(sys2)
>>> sys3 = s22x5[1]
>>> sys3
'Ammonia_dimer_1.0'
>>> atoms3 = create_s22_system(sys3)
>>> get_interaction_energy_s22(sys1)
-0.1375
>>> get_interaction_energy_s22(sys2)
-0.1375
>>> get_interaction_energy_s22(sys3)
-0.1375
>>> get_interaction_energy_s22x5(sys2)
-0.10549743024963291
>>> get_interaction_energy_s22x5(sys3)
-0.1375
>>> get_interaction_energy_s22x5(sys3,correct_offset=False)
-0.1362
>>> get_interaction_energy_s22x5(sys1,dist=1.0)
-0.1375
>>> get_interaction_energy_s22x5(sys1,dist=0.9)
-0.10549743024963291
>>> get_interaction_energy_s22x5(sys1,dist=0.9,correct_offset=False)
-0.1045
>>> get_number_of_dimer_atoms(sys1)
[4, 4]
>>> get_s22x5_distance(sys2)
-0.25040236345454536
>>> get_s22x5_distance(sys3)
0.0

where sys1 is an s22 complex in the original geometry,
while sys2 and sys3 are two different s22x5 geometries
of the exact same complex. It is seen that the interaction
energies for an s22 system and its s22x5 equivalent
(indexed '_1.0') does not neccesarily match
when the energy offset-correction is turned off.
The last two functions are convenience functions,
giving the number of atoms in the two molecules
constituting a dimer and the relative intermolecular
distance in a dimer
(relative to the '1.0' separation, and in Angstrom),
respectively.
