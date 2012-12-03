.. _ase2:

===================================
Porting old ASE-2 code to version 3
===================================

The old Numeric-Python based ASE-2 can coexist with the new
:term:`numpy`-based ASE-3, because the two packages have different
names: *ASE* and *ase* respectively.  Here is an example of combining both::
 
  from ASE.Trajectories.NetCDFTrajectory import NetCDFTrajectory
  traj = NetCDFTrajectory('a.nc')
  loa = traj.GetListOfAtoms(-1)
  # Convert to new style Atoms object:
  from ase import *
  a = Atoms(loa)

The new ASE can actually read old NetCDF trajectory files, so this
would be simpler::

  from ase import *
  a = read('a.nc')

.. note::

   Reading old NetCDF files in the new ASE, works even without having
   the *libnetcdf* and ``Scientific.IO.NetCDF`` libraries installed.



.. index:: ASE2ase

The ASE2ase tool
================

Use the :program:`ASE2ase` tool (source code :svn:`tools/ASE2ase`) to convert old scripts::

  $ ASE2ase oldscript.py
  $ diff -u oldscript.py.bak oldscript.py

Check that the differences look OK.  The conversion tool isn't clever
enough to get everything right, so you will have to do some conversion
manually also.  If you have any problems doing this, then you should
not hesitate to contact the
``campos-devel`` mailing list (see :ref:`mailing_lists`) and ask for help.

