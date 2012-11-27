==================
Atomization energy
==================

The following script will calculate the atomization energy of a
nitrogen molecule:

 .. literalinclude:: N2.py

First, an ``Atoms`` object containing one nitrogen is created and a
fast EMT calculator is attached to it simply as an argument. The total
energy for the isolated atom is then calculated and stored in the
``e_atom`` variable.

The ``molecule`` object is defined, holding the nitrogen molecule at
the experimental bond length. The EMT calculator is then attached to
the molecule and the total energy is extracted into the ``e_molecule``
variable.

Running the script will produce the output::

  Nitrogen atom energy:  4.97 eV
  Nitrogen molecule energy:  0.04 eV
  Atomization energy:  9.90 eV

