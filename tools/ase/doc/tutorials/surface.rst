.. _surface:

================================
Introduction: Nitrogen on copper
================================

This section gives a quick (and incomplete) overview of what ASE can do.

We will calculate the adsorption energy of a nitrogen
molecule on a copper surface.
This is done by calculating the total
energy for the isolated slab and for the isolated molecule. The
adsorbate is then added to the slab and relaxed, and the total energy
for this composite system is calculated. The adsorption energy is
obtained as the sum of the isolated energies minus the energy of the
composite system.

Here is a picture of the system after the relaxation:

.. image:: surface.png

Please have a look at the following script :svn:`doc/tutorials/N2Cu.py`:

.. literalinclude:: N2Cu.py

Assuming you have ASE setup correctly (:ref:`download_and_install`)
run the script::

  python N2Cu.py

Please read below what the script does.

-----
Atoms
-----

The :class:`~ase.atoms.Atoms` object is a collection of atoms.  Here
is how to define a N2 molecule by directly specifying the position of
two nitrogen atoms::

  from ase import Atoms
  d = 1.10
  molecule = Atoms('2N', positions=[(0., 0., 0.), (0., 0., d)])

You can also build crystals using, for example, the lattice module
which returns :class:`~ase.atoms.Atoms` objects corresponding to
common crystal structures. Let us make a Cu (111) surface::

  from ase.lattice.surface import fcc111
  slab = fcc111('Cu', size=(4,4,2), vacuum=10.0)



-----------
Calculators
----------- 

The following :mod:`calculators` can be used with ASE:
:mod:`~calculators.emt`, Asap_, Dacapo_, GPAW_, Siesta_
(Abinit_, Vasp_, and MMTK - work in progress).
  
.. _Asap: http://wiki.fysik.dtu.dk/asap
.. _Dacapo: http://wiki.fysik.dtu.dk/dacapo
.. _GPAW: http://wiki.fysik.dtu.dk/gpaw
.. _Siesta: http://www.icmab.es/siesta
.. _Abinit: http://www.abinit.org
.. _Vasp: http://cms.mpi.univie.ac.at/vasp

In this overview we use the effective medium theory (EMT) calculator,
as it is very fast and hence useful for getting started.

We can attach a calculator to the previously created
:class:`~ase.atoms.Atoms` objects::

  from ase import EMT
  slab.set_calculator(EMT())
  molecule.set_calculator(EMT()) 

and use it to calculate the total energies for the systems by using
the :meth:`~ase.atoms.Atoms.get_potential_energy` method from the
:class:`~ase.atoms.Atoms` class::

  e_slab = slab.get_potential_energy()
  e_N2 = molecule.get_potential_energy()


--------------------
Structure relaxation
--------------------

Let's use the :mod:`QuasiNewton <optimize.qn>` minimizer to optimize the
structure of the N2 molecule adsorbed on the Cu surface. First add the
adsorbate to the Cu slab, for example in the on-top position::
  
  h = 1.85
  add_adsorbate(slab, molecule, h, 'ontop')

In order to speed up the relaxation, let us keep the Cu atoms fixed in
the slab by using :class:`~constraints.FixAtoms` from the
:mod:`~ase.constraints` module. Only the N2 molecule is then allowed
to relax to the equilibrium structure::

  from ase.constraints import FixAtoms
  constraint = FixAtoms(mask=[a.symbol != 'N' for a in slab])
  slab.set_constraint(constraint)

Now attach the :mod:`QuasiNewton <optimize.qn>` minimizer to the
system and save the trajectory file. Run the minimizer with the
convergence criteria that the force on all atoms should be less than
some ``fmax``::

  from ase.optimize import QuasiNewton
  dyn = QuasiNewton(slab, trajectory='N2Cu.traj')
  dyn.run(fmax=0.05)

.. note::

  The general documentation on
  :ref:`structure optimizations <structure_optimizations>` contains
  information about different algorithms, saving the state of an optimizer
  and other functionality which should be considered when performing
  expensive relaxations.

------------
Input-output
------------

Writing the atomic positions to a file is done with the
:func:`~ase.io.write` function::

  from ase.io import write
  write('slab.xyz', slab)

This will write a file in the xyz-format.  Possible formats are:

========  ===========================
format    description
========  ===========================
``xyz``   Simple xyz-format
``cube``  Gaussian cube file
``pdb``   Protein data bank file
``traj``  ASE's own trajectory format
``py``    Python script
========  ===========================

Reading from a file is done like this::

  from ase.io import read
  slab_from_file = read('slab.xyz')

If the file contains several configurations, the default behavior of
the :func:`~ase.io.write` function is to return the last
configuration. However, we can load a specific configuration by
doing::

  read('slab.traj')      # last configuration
  read('slab.traj', -1)  # same as above
  read('slab.traj', 0)   # first configuration


-------------
Visualization
-------------

The simplest way to visualize the atoms is the :func:`~visualize.view`
function::

  from ase.visualize import view
  view(slab)

This will pop up a :mod:`gui` window.  Alternative viewers can be used
by specifying the optional keyword ``viewer=...`` - use one of
'ase.gui', 'gopenmol', 'vmd', or 'rasmol'. (Note that these alternative
viewers are not a part of ASE and will need to be installed by the user
separately.) The VMD viewer can take an optional ``data`` argument to
show 3D data::

  view(slab, viewer='VMD', data=array)


------------------
Molecular dynamics
------------------

Let us look at the nitrogen molecule as an example of molecular
dynamics with the :class:`VelocityVerlet <md.verlet.VelocityVerlet>`
algorithm. We first create the :class:`VelocityVerlet
<md.verlet.VelocityVerlet>` object giving it the molecule and the time
step for the integration of Newton's law. We then perform the dynamics
by calling its :meth:`run` method and giving it the number of steps to
take::

  from ase.md.verlet import VelocityVerlet
  from ase import units
  dyn = VelocityVerlet(molecule, dt=1.0 * units.fs)
  for i in range(10):
     pot = molecule.get_potential_energy()
     kin = molecule.get_kinetic_energy()
     print '%2d: %.5f eV, %.5f eV, %.5f eV' % (i, pot + kin, pot, kin)
     dyn.run(steps=20)


