.. module:: calculate

=========
Calculate
=========

Set calculator
--------------

Allows ag to choose a calculator for internal computations (see
below). Different density functional codes and force fields, as well
as the EMT calculator are available. For the FHI-aims and VASP
calculators, it is also possible to export an entire set of input
files. 

Energy and forces
-----------------

Invokes the currently set calculator and provides energies and
optional forces for all atoms. 

Energy minimization
-------------------

Runs an ASE relaxation using the currently selected calculator with a
choice of relaxation algorithm and convergence criteria. Great for
quickly (pre-)relaxing a molecule before placing it into a bigger
system. 

Scale system
------------

Stub
