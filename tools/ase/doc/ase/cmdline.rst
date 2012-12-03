.. highlight:: bash

.. index:: ase

.. _command line tools:

==================
Command line tools
==================

The :program:`ase` program can be used to do calculations with
:ref:`ASE supported calculators<supported calculators>` on the command
line without having to write a Python script.  The syntax is::

    $ ase [calculator] [task] [options] system(s)

The name of the calculator must be lower case and will default to
:mod:`EMT <emt>`.  The task must be ``molecule`` or ``bulk``.  There are
several ways to specify the system or systems to perform the
calculations on:

* Chemical names: ``H2O`` or ``Fe``.
  Default :ref:`molecule<molecules-section>` definitions are used.
* Range of chemical symbols: ``Sc-Zn`` (3d-metals)
* Special names: ``G2``, ``G2-1`` or ``S22``
* File names: ``benzene.xyz`` or ``slab.traj``

The exact meaning of these names will depend on the task.

Simple examples::

    $ ase emt H2 --relax=0.01
    $ ase abinit bulk Si -a 5.5 -p ecut=150 -k 4,4,4


Command line options
====================

General options:

-h, --help          Show help message and exit.
-t TAG, --tag=TAG   String tag added to filenames.
-M <M1,M2,...>, --magnetic-moment=<M1,M2,...>
                    Magnetic moment(s).  Use "-M 1" or "-M 2.3,-2.3".
-G, --gui           Pop up ASE's GUI.
-s, --write-summary
                    Write summary.
--slice=<start:stop:step>
                    Select subset of calculations using Python slice
                    syntax.  Use "::2" to do every second calculation and
                    ":-5" to do the last five.
-w FILENAME, --write-to-file=FILENAME
                    Write configuration to file.
-i, --interactive-python-session
                    Run calculation inside interactive Python session.  A
                    possible $PYTHONSTARTUP script will be imported and
                    the "atoms" variable refers to the Atoms object.
-l, --use-lock-files
                    Skip calculations where the json lock-file or result file
                    already exists.
-R FMAX, --relax=FMAX
                    Relax internal coordinates using L-BFGS algorithm.
-F <N,x>, --fit=<N,x>
                    Find optimal bondlength and vibration frequency for
                    dimer molecules or optimal volume and bulk modulus for
                    bulk systems using N points and a variation from -x %
		    to +x % for the bondlength or lattice constants.
--constrain-tags=<T1,T2,...>
                    Constrain atoms with tags T1, T2, ...
-k <K1,K2,K3>, --monkhorst-pack=<K1,K2,K3>
                    Monkhorst-Pack sampling of BZ.  Example: "4,4,4":
                    4x4x4 k-points, "4,4,4g": same set of k-points shifted
                    to include the Gamma point.
--k-point-density=K_POINT_DENSITY
                    Density of k-points in Å.
-p <key=value,...>, --parameters=<key=value,...>
                    Comma-separated key=value pairs of calculator specific
                    parameters.

Options specific to the molecule task:

-v VACUUM, --vacuum=VACUUM
                    Amount of vacuum to add around isolated systems (in
                    Angstrom).
--unit-cell=CELL
                    Unit cell.  Examples: "10.0" or "9,10,11" (in
                    Angstrom).
--bond-length=BOND_LENGTH
                    Bond length of dimer in Angstrom.
--atomize           Calculate Atomization energies.

Options specific to the bulk task:

-x CRYSTAL_STRUCTURE, --crystal-structure=CRYSTAL_STRUCTURE
                    Crystal structure.
-a LATTICE_CONSTANT, --lattice-constant=LATTICE_CONSTANT
                    Lattice constant in Å.
--c-over-a=C_OVER_A
                    c/a ratio.
-O, --orthorhombic  Use orthorhombic unit cell.
-C, --cubic         Use cubic unit cell.
-r REPEAT, --repeat=REPEAT
                    Repeat unit cell.  Use "-r 2" or "-r 2,3,1".


Molecules
=========

Example::

    $ ase abinit H2 -p ecut=200,xc=LDA -F 5,1 --atomize

This will calculate the energy of a :mol:`H_2` molecule using
:mod:`Abinit <abinit>` with a planewave cutoff of 200 eV and the LDA
XC-functional.  A fit using 5 points and a variation of the bond
length from -1 % to +1 % is made and in addition the energy of a
single hydrogen atom is also calculated.

Results are written to json files and can be analysed with::

    $ ase abinit H H2 -s
        name          E       E-E0         d0        hnu         Ea        Ea0
                     eV         eV        Ang        meV         eV         eV
          H2    -29.703      0.022      0.770    556.096      4.852      4.873
           H    -12.426  


.. note::

    The json files are simple text files that can be more or less'ed
    or pretty printed with ``python -m json.tool
    H2-molecule-abinit.jon``.


Bulk systems
============

Example::

    $ ase bulk Ni Cu Pd Ag Pt Au -F 5,1

Here we used the default EMT potential and the result is::

    $ ase bulk Ni Cu Pd Ag Pt Au -s
        name          E       E-E0         V0          B
                     eV         eV      Ang^3        GPa
          Ni     -0.009      0.005     10.600    175.978
          Ag      0.002      0.002     16.775    100.151
          Pt     -0.000      0.000     15.080    278.087
          Au      0.003      0.003     16.684    173.868
          Pd      0.000      0.001     14.588    179.105
          Cu     -0.006      0.001     11.565    134.439


More examples
-------------

Anti-ferromagnetic bcc iron::

    $ ase vasp bulk -x bcc Fe -C -M 2.3,-2.3 -p xc=PBE -k 8,8,8

Bulk silicon (128 atoms, `\Gamma` point only)::

    $ ase abinit bulk Si -r 4,4,4 -k 1,1,1 -a 5.46

Bulk aluminum in orthorhombic cell with LDA and fixed rmt::

    $ ase elk bulk --orthorhombic Al -k 4,4,4 -a 4.05 -p "swidth=0.1,rmt={'Al': 1.9}"

Batch jobs
==========

Suppose you want to run a large number of similar calculations like
relaxing the structure of all the molecules in the :ref:`G2-1 database
<molecular data>`.  You could do that by submitting this job to your
compute cluster::

    $ ase gpaw G2-1 -v 6.0 -p xc=vdW-DF,h=0.18 -R 0.02

The molecule task will expand ``G2-1`` to a lot of molecules, so it
makes sense to use :option:`-l` option (:option:`--use-lock-files`)
and submit the same job many times. A lock file will be created
for each started calculation and calculations with existing lock file skipped. 
Moreover the calculations can be run in parallel
(if parallel version of GPAW is installed)::

    $ mpiexec gpaw-python `which ase` gpaw G2-1 -v 6.0 -p xc=vdW-DF,h=0.18 -R 0.02 -l


Making your own tasks
=====================

FCC clusters with 13 atoms
--------------------------

Put this in :file:`m13.py`:

.. literalinclude:: m13.py
    :language: python

Then do this::

    $ ase m13.py Pt -R 0.01

The relaxed EMT bondlength of 2.62 Å can be extracted from the
created trajectory like this::

    $ ag -tg "e,fmax,d(0,9)" Pt-m13-emt.traj
    10.9824302706 2.27575724833 2.72
    10.0805658256 1.45353946744 2.68
    9.61654400875 0.447140179352 2.64
    9.57510700742 0.0656434401881 2.6222281202
    9.57424531861 0.00239771758341 2.62450316817

or like this:

.. highlight:: python

>>> from ase.io import read
>>> pt = read('Pt-m13-emt.traj')
>>> pt.get_distance(0, 9)
2.6245031681662452


Convergence test
----------------

See `convergence.py`_.

.. _convergence.py: https://trac.fysik.dtu.dk/projects/gpaw/browser/trunk/gpaw/
                    tasks/convergence.py


To be done
==========

* Optimize c/a ratio.
* Implement different way of cell sampling for eos:
  [a0 + s for s in [x * np.array(range(- N/2 + 1, N/2 + 1))]]
  where x is sampling step length, N number of steps.
  Current way of sampling gives different length of sampling interval
  depending on the lattice constant guess a0.
* Write results to file (pickel, csv, cmr (db), traj, ...) per system,
  together with json file!
* Split off EnergyTask from Task.
* Set correct magnetic moments for atoms. DONE
* Allow setting charges in ase.tasks.task
* Check occupation numbers and requested total magnetic moments
  for molecules/atoms. DONE
* Add --exclude option.
* Relax cell.
* Optimize first then fit.
* Behavior of -w option?
* Reaction task?
* Rethink analysis and summary stuff:
   * it would be nice to calculate for example cohesive energies on 
     the command line, i.e. using species belonging to different tasks
   * analysis should be more modular: one may want for example to calculate
     zpe energies for adsorption systems including molecules and surfaces
     and print the zpe correction for a given reaction.
* ase.tasks.main - print complete architecture string in case of error
  (like in ase/test/__init__.py)

