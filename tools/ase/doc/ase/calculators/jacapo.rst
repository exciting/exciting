.. module::  jacapo
   :synopsis: ASE python interface for Dacapo

==========================================================
Jacapo - ASE python interface for Dacapo
==========================================================

Introduction
============

Jacapo_ is an ASE interface for Dacapo_ and fully compatible with ASE. It 
replaces the old Dacapo interface using Numeric python and ASE2.
The code was originally developed by John Kitchin and detailed documentation
as well as many examples are available online:

http://gilgamesh.cheme.cmu.edu/doc/software/jacapo/index.html

Jacapo is included as an optional calculator in ASE and small differences to the
above documentation may occur.

.. _Jacapo: http://gilgamesh.cheme.cmu.edu/doc/software/jacapo/index.html
.. _Dacapo: http://wiki.fysik.dtu.dk/dacapo

Jacapo Calculator
================= 

The Jacapo calculator is automatically installed with ase and can be imported using::

  from ase.calculators.jacapo import *

.. class:: Jacapo()
    
Here is a list of available keywords to initialize the calculator:

============== ============ =====================================
keyword        type         description
============== ============ =====================================
``nc``         ``str``      Output NetCDF file, or input file if nc already exists.
``outnc``      ``str``      Output file. By default equal to nc.
``atoms``      ``object``   ase atoms object
``pw``         ``float``    Planewave cutoff in eV
``dw``         ``float``    Density cutoff in eV
``xc``         ``str``      Exchange-correlation functional. One of ['PZ','VWN','PW91','PBE','RPBE','revPBE']
``nbands``     ``int``      Number of bands
``ft``         ``float``    Fermi temperature
``kpts``       ``list``     K-point grid, e.g. kpts = (2,2,1)
``spinpol``    ``boolean``  Turn on/off spin-polarization
``fixmagmom``  ``str``      Magnetic moment of the unit cell
``symmetry``   ``boolean``  Turn on/off symmetry reduction
``stress``     ``boolean``  Turn on/off stress calculation
``dipole``     ``boolean``  Turn on/off dipole correction
``stay_alive`` ``boolean``  Turn on/off stay alive
``debug``      ``int``      Set debug level (0=off, 10=extreme)
============== ============ =====================================

Example
=======

Here is an example of how to calculate the total energy of a CO molecule::
        
  #!/usr/bin/env python
  from ase import *
  from ase.calculators.jacapo import *

  CO = data.molecules.molecule('CO')
  CO.set_cell([6,6,6])
  CO.center()

  calc = Jacapo(nc='CO.nc',
                atoms=CO,
                pw=300,
                nbands=8)
  
  print CO.get_potential_energy()
  

Restarting from an old Calculation
==================================

With Jacapo it is very easy to restart an old calculation. All necessary information
(including the constraints) is stored in the out.nc file. To continue for example a
geometry optimization do something like this::

  calc = Jacapo('old.nc',stay_alive=True)
  atoms = calc.get_atoms()
  dyn = QuasiNewton(atoms,logfile='qn.log')
  dyn.run(fmax=0.05)

Note, that the stay_alive flag is not stored in the .nc file and must be set when the
calculator instance is created.


