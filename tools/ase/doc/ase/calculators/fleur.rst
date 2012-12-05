.. module:: fleur

=====
FLEUR
=====

Introduction
============

FLEUR_ is a density-functional theory code which uses the full
potential linearized augmented plane-wave (FLAPW) method. FLEUR_ can
be applied to any element in the periodic table, and the code is well
suited especially to surfaces and magnetic materials.


.. _FLEUR: http://www.flapw.de



Environment variables
=====================

In order to use FLEUR through ASE, two environment variables have to
be set. :envvar:`FLEUR_INPGEN` should point to the simple input
generator of FLEUR, and :envvar:`FLEUR` to the actual executable. Note
that FLEUR has different executables e.g. for cases with and without
inversion symmetry, so the environment variable has to be set accordingly.
As an example, the variables could be set like:

.. highlight:: bash
 
::

  $ export FLEUR_INPGEN=$HOME/fleur/v25/inpgen.x
  $ export FLEUR=$HOME/fleur/v25/fleur.x

.. highlight:: python

or:

.. highlight:: bash
 
::

  $ export FLEUR="mpirun -np 32 $HOME/fleur/v25/fleur.x"

.. highlight:: python

for parallel calculations.

FLEUR Calculator
================

Currently, limited number of FLEUR parameters can be set via the ASE interface 
Below follows a list of supported parameters

===============  =========  ==============  ============================
keyword          type       default value   description
===============  =========  ==============  ============================
``xc``           ``str``    ``'LDA'``       XC-functional. Must be one 
                                            of 'LDA', 'PBE', 'RPBE'
``kpts``         *seq*      `\Gamma`-point  **k**-point sampling
``convergence``  *dict*     ``{'energy':    Convergence criteria (meV)
                            0.0001}``
``width``        ``float``                  Width of Fermi smearing (eV)
``kmax``         ``float``                  Plane-wave cut-off (a.u.)
``mixer``        *dict*                     Mixing parameters 'imix',
                                            'alpha', and 'spinf'
``maxiter``      ``int``    40              Maximum number of SCF steps
``maxrelax``     ``int``    20              Maximum number of relaxation 
                                            steps
``workdir``      ``str``    Current dir     Working directory for the
                                            calculation
===============  =========  ==============  ============================

*seq*: A sequence of three ``int``'s.
*dict*: A dictionary


Spin-polarized calculation
==========================

If the atoms object has non-zero magnetic moments, a spin-polarized calculation
will be performed by default.

Utility functions
=================

As only a subset of FLEUR parameters can currently be specified
through ASE interface, the interface defines some utility functions
for cases where manual editing of the FLEUR input file ``inp`` is
necessary.

.. automethod:: ase.calculators.fleur.FLEUR.write_inp
.. automethod:: ase.calculators.fleur.FLEUR.initialize_density
.. automethod:: ase.calculators.fleur.FLEUR.calculate
.. automethod:: ase.calculators.fleur.FLEUR.relax

Examples
========

Lattice constant of fcc Ni

.. literalinclude:: fcc_Ni_fleur.py

See the :ref:`equation of states tutorial <eos>` for analysis of the results.

