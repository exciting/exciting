.. module:: abinit

======
ABINIT
======

Introduction
============

ABINIT_ is a density-functional theory code
based on pseudopotentials and a planewave basis.


.. _ABINIT: http://www.abinit.org



Environment variables
=====================

**Note**: setting environment variables is necessary
only if you configure your ABINIT/ASE environment from scratch.

You need to write a script called :file:`abinit.py` containing
something like this::

  import os
  abinit = '/usr/bin/abinis'
  exitcode = os.system('%s < %s.files > %s.log' % (abinit, label, label))

The environment variable :envvar:`ABINIT_SCRIPT` must point to that file.

A directory containing the pseudopotential files :file:`.fhi` is also
needed, and it is to be put in the environment variable
:envvar:`ABINIT_PP_PATH`.

Set both environment variables in your in your shell configuration file:

.. highlight:: bash
 
::

  $ export ABINIT_SCRIPT=$HOME/bin/abinit.py
  $ export ABINIT_PP_PATH=$HOME/mypps

.. highlight:: python



ABINIT Calculator
================= 

The default parameters are very close to those that the ABINIT Fortran
code uses.  These are the exceptions:

.. class:: Abinit(label='abinit', xc='LDA', pulay=5, mix=0.1)
    
Here is a detailed list of all the keywords for the calculator:

============== ========= ================  =====================================
keyword        type      default value     description
============== ========= ================  =====================================
``kpts``       ``list``  ``[1,1,1]``       Monkhorst-Pack k-point sampling
``nbands``     ``int``   ``None``          Number of band. May be omitted.
``nstep``      ``int``   ``None``          Number of self-consistent field STEPS.
``ecut``       ``float`` ``None``          Planewave cutoff energy in eV (default: None)
``xc``         ``str``   ``'LDA'``         Exchange-correlation functional.
``pulay``      ``int``   ``5``             Number of old densities to use for
                                           Pulay mixing
``mix``        ``float`` ``0.1``           Pulay mixing weight 
``width``      ``float`` ``0.04 Ha``       Fermi-distribution width in eV (default: 0.04 Ha)
``charge``     ``float`` ``0``             Total charge of the system (default: 0)
``label``      ``str``   ``'abinit'``      Name of the output file
``toldfe``     ``float`` ``1.0e-6``        TOLerance on the DiFference of total Energy
============== ========= ================  =====================================

A value of ``None`` means that ABINIT's default value is used.

**Warning**: abinit does not specify a default value for
``Planewave cutoff energy in eV`` - you need to set them as in the example at thei bottom of the page, otherwise calculation will fail.
Calculations wihout k-points are not parallelized by default
and will fail! To enable band paralellization specify ``Number of BanDs in a BLOCK`` 
(``nbdblock``) as `Extra parameters`_ -
see `<http://www.abinit.org/Infos_v5.2/tutorial/lesson_parallelism.html>`_.

Extra parameters
================

The ABINIT code reads the input parameters for any calculation from a 
:file:`.in` file and :file:`.files` file.
This means that you can set parameters by manually setting 
entries in this input :file:`.in` file. This is done by the syntax:

>>> calc.set_inp('name_of_the_entry', value)

For example, the ``ndtset`` can be set using

>>> calc.set_inp('ndtset', 2)

The complete list of keywords can be found in the official `ABINIT
manual`_.

.. _ABINIT manual: http://www.abinit.org/Infos_v5.4/input_variables/keyhr.html



Pseudopotentials
================

Pseudopotentials in the ABINIT format are available on the
`pseudopotentials`_ website.
A database of user contributed pseudopotentials is also available there.

.. _pseudopotentials: http://www.abinit.org/Psps/?text=psps



Example 1
=========

Here is an example of how to calculate the total energy for bulk Silicon:
        
.. literalinclude:: Si_abinit.py
