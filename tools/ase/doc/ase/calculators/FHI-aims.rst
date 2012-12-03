.. module:: FHI-aims

========
FHI-aims
========

Introduction
============

FHI-aims_ is a all-electron full-potential density functional theory 
code using a numeric local orbital basis set. This interface provides
all that should be required to run FHI-aims_ from within ASE.

.. _FHI-aims: http://www.fhi-berlin.mpg.de/aims/

Running the Calculator
====================== 

The default initialization command for the FHI-aims calculator is 

.. class:: Aims(output_template = 'aims', track_output = False)

In order to run a calculation, you have to ensure that at least the 
following ``str`` variables are specified, either in the initialization 
or as shell variables:

===============  ====================================================
keyword          description
===============  ====================================================
``run_command``   The full command required to run FHI-aims from 
		  a shell, including anything to do with an MPI
		  wrapper script and the number of tasks.
		  An alternative way to set this command is via the 
		  shell variable ``AIMS_COMMAND``, which is checked
		  upon initialization and when starting a run. 
``species_dir``   Directory where the species defaults are located 
		  that should be used. Can also be specified with 
		  the system variable ``AIMS_SPECIES_DIR``.
``xc``            The minimal physical specification: what kind of 
		  calculation should be done. 
===============  ====================================================

In addition, you might want to specify at least one of self-consistency 
accuracy commands (see below) in order to avoid an excessively long 
calculation. 

Two general options might come in useful to postprocess the output:

===================  ====================================================
keyword              description
===================  ====================================================
``output_template``  Base name for the output, in case the calculator
		     is called multiple times within a single script. 
``track_output``     ``True/False`` - if ``True`` all the output files
		     will be kept, while the number of calls to the 
		     calculator is encoded in the output file name. 
===================  ====================================================

List of keywords
================

This is a non-exclusive list of keywords for the ``control.in`` file 
that can be addresses from within ASE. The meaning for these keywords is 
exactly the same as in FHI-aims, please refer to its manual for help on 
their use. 

One thing that should be mentioned is that keywords with more than
one option have been implemented as lists, eg. 
``k_grid=(12,12,12)`` or ``relativistic=('atomic_zora','scalar')``. 
In those cases, specifying a single string containing all the options is also possible. 

None of the keywords have any default within ASE,but do check the defaults
set by FHI-aims. If there is a keyword that you would 
like to set and that is not yet implemented here, it is trivial to add 
to the first few lines of the aims calculator in the file 
ASE/ase/calculators/aims.py .

Describing the basic physics of the system:

============================  ======
keyword                       type 
============================  ======
``xc``			      str   
``charge``                    float
``spin``		      str
``relativistic``	      list
``use_dipole_correction``     bool
``vdw_correction_hirshfeld``  str
``k_grid``		      list
============================  ======

Driving relaxations and molecular dynamics:

============================  ======
keyword                       type 
============================  ======
``relax_geometry``	      list
``max_relaxation_steps``      int
``n_max_pulay``  	      int
``sc_iter_limit``	      int
``restart_relaxations``	      bool
``MD_run``		      list
``MD_schedule``		      list
``MD_segment``		      list
============================  ======

Output options:

============================  ========
keyword                       type 
============================  ========
``output_level``	      str
``output``		      list
``cubes``                     AimsCube
============================  ========

See below for a description of the volumetric cube file output
interface AimsCube

Keywords for accuracy settings:

============================  ======
keyword                       type 
============================  ======
``sc_accuracy_eev``	      exp
``sc_accuracy_etot``	      exp
``sc_accuracy_forces``	      exp
``sc_accuracy_rho``	      exp
``compute_forces``	      bool
============================  ======

Keywords to adjust the SCF-cycle

============================  ======
keyword                       type 
============================  ======
``charge_mix_param``	      float	
``prec_mix_param``	      float
``spin_mix_param``	      float
``KS_method``		      str
``restart``		      str
``restart_read_only``	      str
``restart_write_only``	      srt
``preconditioner``	      list
``mixer``		      str
``empty_states``	      int	
``ini_linear_mixing``	      int
``mixer_threshold``	      list
``occupation_type``	      list
============================  ======

Note:: 
   
 Any argument can be changed after the initial construction of the
 calculator, simply by setting it with the method 

   >>> calc.set( keyword=value )

Volumetric Data Output
======================

The class

.. class:: AimsCube(origin=(0,0,0),edges=[(0.1,0.0,0.0),(0.0,0.1,0.0),(0.0,0.0,0.1)],points=(50,50,50),plots=None)

describes an object that takes care of the volumetric
output requests within FHI-aims. An object of this type can 
be attached to the main Aims() object as an option. 

The possible arguments for AimsCube are:

============================  ========
keyword                       type 
============================  ========
``origin``		      list
``edges``		      3x3-array
``points``		      list
``plots``		      list
============================  ========

The possible values for the entry of plots 
are discussed in detail in the FHI-aims manual, 
see below for an example.

Example
=======

As an example, here is a possible setup to obtain 
the geometry of a water molecule:

.. literalinclude:: H2O_aims.py
