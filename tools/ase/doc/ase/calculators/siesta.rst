.. module:: siesta

======
SIESTA
======

Introduction
============

SIESTA_ is a density-functional theory code for very large systems
based on atomic orbital (LCAO) basis sets.


.. _SIESTA: http://www.uam.es/siesta/



Environment variables
=====================

You need to write a script called :file:`run_siesta.py` containing
something like this::

  import os
  siesta = 'siesta_2.0'
  exitcode = os.system('%s < %s.fdf > %s.txt' % (siesta, label, label))

The environment variable :envvar:`SIESTA_SCRIPT` must point to that file.

A directory containing the pseudopotential files :file:`.vps` is also
needed, and it is to be put in the environment variable
:envvar:`SIESTA_PP_PATH`.

Set both environment variables in your shell configuration file:

.. highlight:: bash
 
::

  $ export SIESTA_SCRIPT=$HOME/siesta/run_siesta.py
  $ export SIESTA_PP_PATH=$HOME/mypps

.. highlight:: python



SIESTA Calculator
================= 

The default parameters are very close to those that the SIESTA Fortran
code uses.  These are the exceptions:

.. class:: Siesta(label='siesta', xc='LDA', pulay=5, mix=0.1)
    
Here is a detailed list of all the keywords for the calculator:

============== ========= ============= =====================================
keyword        type      default value description
============== ========= ============= =====================================
``kpts``       ``list``  ``[1,1,1]``   Monkhorst-Pack k-point sampling
``nbands``     ``int``   ``None``      Number of bands 
``meshcutoff`` ``float`` ``None``      Mesh cut-off energy in eV 
``basis``      ``str``   ``None``      Type of basis set ('sz', 'dz', 'szp',
                                       'dzp') 
``xc``         ``str``   ``'LDA'``     Exchange-correlation functional.
``pulay``      ``int``   ``5``         Number of old densities to use for
                                       Pulay mixing
``mix``        ``float`` ``0.1``       Pulay mixing weight 
``width``      ``float`` ``None``      Fermi-distribution width in eV
``charge``     ``float`` ``None``      Total charge of the system
``label``      ``str``   ``'siesta'``  Name of the output file
============== ========= ============= =====================================

A value of ``None`` means that SIESTA's default value is used.



Extra FDF parameters
====================

The SIESTA code reads the input parameters for any calculation from a 
:file:`.fdf` file. This means that you can set parameters by manually setting 
entries in this input :file:`.fdf` file. This is done by the syntax:

>>> calc.set_fdf('name_of_the_entry', value)

For example, the ``EnergyShift`` can be set using

>>> calc.set_fdf('PAO.EnergyShift', 0.01 * Rydberg)

The complete list of the FDF entries can be found in the official `SIESTA
manual`_.

.. _SIESTA manual: http://www.uam.es/departamentos/ciencias/fismateriac/siesta



Customized basis-set
====================

Standard basis sets can be set by the keyword ``basis`` directly
in the Siesta calculator object. If a customized basis is needed, it 
can be set as an FDF entry, as explained in the previous section.

As an example, we generate a triple-zeta triple-polarized (TZTP)
basis for Au. Since the valence states are 6s and 5d, we will have
3 zeta orbitals for l=0 and 3 for l=2 plus 3 polarization orbitals
for l=1. The basis can be defined by

>>> value = ["""Au   2   split  0.00  #label, num. of l-shells,type,charge
>>>         0   3   P    3            #l,nzeta,'P'(opt):pol.functions,npolzeta
>>>         0.00   0.00   0.00        #rc of basis functions for each zeta function
>>>                                   #0.00  => rc determined by PAO.EnergyShift
>>>         2   3                     #l,nzeta
>>>         0.00   0.00   0.00"""]    #rc

>>> calc.set_fdf('PAO.Basis',value=value)

The default basis set generation fails for Pt for some versions of
Siesta. If this happens, you should specify the basis set
manually. This is exemplified below.

For Pt, using ``calc.set_fdf('PAO.EnergyShift', 0.1 * eV)`` is usually
resonable, and a single-zeta polarized basis set can be specified by
writing::

  # Define SZP basis set for Pt
  calc.set_fdf('PAO.Basis',
               ["""\
  Pt   2         # Species label, number of l-shells
  n=6  0  1 P    # n, l, Nzeta, Polarization, NzetaPol
  0.             # 0.0 => default [6.982 \n 1.000]
  n=5  2  1      # n, l, zeta
  0."""])        # 0.0 => default [5.172 \n 1.000]

A double-zeta polarized basis set for Pt may be specified by::

  # Define DZP basis set for Pt
  calc.set_fdf('PAO.Basis',
               ["""\
  Pt 2 split 0.00  # Species label, number of l-shells
  n=6 0 2 P 1      # n, l, Nzeta, Polarization, NzetaPol
  0.00 0.00        # 0.0 => default [6.982  5.935 \n 1.000  1.000]
  n=5 2 2          # n, l, zeta
  0.00 0.00"""])   # 0.0 => default [5.172  3.060 \n 1.000  1.000]

You can also reuse the basis set of a previous calculation, by copying
the .ion files to the new location, and set the ``User.Basis`` tag to
``True``::

  # Load basis from previous calc (*.ion files)
  calc.set_fdf('User.Basis', True)

Warning: Specifying a basis set manually in any way will, for some
obscure reason, make Siesta crash if you have ghost atoms!



Pseudopotentials
================

Pseudopotential files in the ``.psf`` or ``.vps`` formats are needed. 
Pseudopotentials generated from the ABINIT code and converted to 
the SIESTA format are available in the `SIESTA`_ website . A database of user 
contributed pseudopotentials is also available there.

You can also find an on-line pseudopotential generator_ from the
OCTOPUS code.

.. _generator: http://www.tddft.org/programs/octopus/wiki/index.php/Pseudopotentials



Example
=======

Here is an example of how to calculate the total energy for bulk Silicon,
using a double-zeta basis generated by specifying a given energy-shift::
        
  #!/usr/bin/env python
  from ase import *
  
  a0 = 5.43
  bulk = Atoms('Si2', [(0, 0, 0),
                       (0.25, 0.25, 0.25)],
               pbc=True)
  b = a0 / 2
  bulk.set_cell([(0, b, b),
                 (b, 0, b),
                 (b, b, 0)], scale_atoms=True)
  
  calc = Siesta(label='Si',
                xc='PBE',
                meshcutoff=200 * Ry,
                basis='dz',
                mix=0.01,
                kpts=[10, 10, 10])
   
  calc.set_fdf('PAO.EnergyShift', 0.01 * Ry)
  bulk.set_calculator(calc)
  e = bulk.get_potential_energy()

Here, the only input information on the basis set is, that it should
be double-zeta (``basis='dz'``) and that the confinement potential
should result in an energy shift of 0.01 Rydberg (the
``PAO.EnergyShift`` fdf tag). Sometimes it can be necessary to specify
more information on the basis set. For example, the default basis set
generation fails for Pt for some versions of Siesta. To fix this, you
*must* specify the basis set manually. Manual basis set specifications
are described in `Customized basis-set`_.



Restarting from an old Calculation
==================================

If you want to rerun an old SIESTA calculation, made using the ASE
interface or not, you can set the fdf tag ``UseSaveData`` to
``True``. This is equivalent to setting both ``DM.UseSaveDM`` and
``MD.UseSaveXV`` to True, i.e. it will reuse the the density matrix,
and the atomic coordinates (and unit cell) of the previous
calculation.  Note that the Siesta jobname (the ``label`` keyword in
the ASE interface) must be identical to the jobname of the old
calculation.
