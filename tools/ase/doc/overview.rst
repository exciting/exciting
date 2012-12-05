.. _overview:

========
Overview
========

ASE is an Atomistic Simulation Environment written in the
Python_ programming language with the aim of setting up, stearing, and
analyzing atomistic simulations. The ASE has been constructed with a
number of "design goals" that make it:


- **Easy to use**:

  Setting up an atomistic total energy calculation or molecular
  dynamics simulation with ASE is simple and straightforward.  ASE can
  be used via a :mod:`graphical user interface <gui>`, :ref:`command
  line tools` and the Python language.  Python scripts are
  easy to follow (see :ref:`python_info` for a short introduction).
  It is simple for new users to get access to all of the functionality
  of ASE.

- **Flexible**:

  Since ASE is based on the Python scripting language it is possible
  to perform very complicated simulation tasks without any code modifications.
  For example, a sequence of calculations may be performed with
  the use of simple "for-loop" constructions. There exist ASE modules for 
  performing many standard simulation tasks.

- **Customizable**:

  The Python code in ASE is structured in modules intended for
  different purposes. There are :mod:`calculators` for calculating
  energies, forces and stresses, :mod:`md` and :mod:`optimize` modules
  for controlling the motion of atoms, :mod:`constraint <constraints>`
  objects and filters for performing :mod:`nudged-elastic-band <neb>`
  calculations etc. The modularity of the object-oriented code make it 
  simple to contribute new functionality to ASE.

- **Pythonic**:

  It fits nicely into the rest of the Python world with
  use of the popular :term:`NumPy` package for numerical work
  (see :ref:`numpy` for a short introduction). The
  use of the Python language allows ASE to be used both interactively
  as well as in scripts.

- **Open to participation**:

  The CAMPOS Atomic Simulation Environment is released under the GNU
  Lesser General Public License version 2.1 or any later version.  See
  the files :trac:`COPYING` and :trac:`COPYING.LESSER` which accompany
  the downloaded files, or see the license at GNU's web server at
  http://www.gnu.org/licenses/.  Everybody is invited to
  participate in using and :ref:`developing the code <devel>`.

.. _Python: http://www.python.org
