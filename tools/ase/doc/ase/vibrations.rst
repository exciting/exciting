.. module:: vibrations

Vibration analysis
------------------

You can calculate the vibrational modes of a an
:class:`~ase.atoms.Atoms` object in the harmonic approximation using
the :class:`~ase.vibration.Vibrations`.

.. autoclass:: ase.vibrations.Vibrations
   :members:

name is a string that is prefixed to the names of all the files
created. atoms is an Atoms object that is either at a
fully relaxed ground state or at a saddle point. freeatoms is a
list of atom indices for which the vibrational modes will be calculated,
the rest of the atoms are considered frozen. displacements is a
list of displacements, one for each free atom that are used in the
finite difference method to calculate the Hessian matrix. method is -1
for backward differences, 0 for centered differences, and 1 for
forward differences.

.. warning::
   Using the `dacapo` caculator you must make sure that the symmetry
   program in dacapo finds the same number of symmetries for the
   displaced configurations in the vibrational modules as found in
   the ground state used as input.
   This is because the wavefunction is reused from one displacement
   to the next.
   One way to ensure this is to tell dacapo not to use symmetries.

   This will show op as a python error 'Frames are not aligned'.
   This could be the case for other calculators as well.


You can get a NetCDF trajectory corresponding to a specific mode by
using:

>>> mode=0
>>> vib.create_mode_trajectory(mode=mode,scaling=5)

This will create a NetCDF trajectory file `CO_vib_mode_0.traj`,
corresponding to the highest frequency mode.
`scaling` is an option argument, that will give the amplitude of
the mode, default is 10.

