.. _structure_optimizations:

======================
Structure optimization
======================

The optimization algorithms can be roughly devided into local
optimization algorithms which find the next local minimum and
global optimization algorithms that try to find the global
minimum (a much harder task).


.. seealso::

    `Performance test
    <https://wiki.fysik.dtu.dk/gpaw/devel/ase_optimize/ase_optimize.html>`_ for all
    ASE optimizers.
 


Local optimization
==================
.. module:: optimize
   :synopsis: Structure Optimization

There are currently 5 different optimization algorithms available:
``BFGS``, ``LBFGS``, ``BFGSLineSearch``, ``LBFGSLineSearch``,
``MDMin``, and ``FIRE``.

``MDMin`` and ``FIRE`` both use Newtonian dynamics with added
friction, to converge to an energy minimum, whereas the first 3 are of
the quasi-Newton type, where the forces of consecutive steps are used
to dynamically update a Hessian describing the curvature of the
potential energy landscape.  You can use the ``QuasiNewton`` synonym
for ``BFGSLineSearch`` because this algorithm is in many cases the optimal one
of the three quasi-Newton algorithms.

All optimizer classes have the following structure::

  class Optimizer:
      def __init__(self, atoms, restart=None, logfile=None):
      def run(self, fmax=0.05, steps=100000000):
      def get_number_of_steps():

The convergence criterion is that the force on all individual atoms
should be less than *fmax*:

.. math:: \max_a |\vec{F_a}| < f_\text{max}


BFGS
----
.. module:: optimize.qn
   :synopsis: Quasi-Newton

The ``BFGS`` object is one of the minimizers in the ASE
package.  Let's try to use it to optimize the structure of a water
molecule.  We start with the experimental geometry::

  from ase import *
  import numpy as np
  d = 0.9575
  t = pi / 180 * 104.51
  water = Atoms('H2O',
                positions=[(d, 0, 0),
                           (d * np.cos(t), d * np.sin(t), 0),
                           (0, 0, 0)],
                calculator=EMT())
  dyn = BFGS(water)
  dyn.run(fmax=0.05)
  BFGS:   0  16:14:26        6.445801      51.6847
  BFGS:   1  16:14:26        2.418583      27.2946
  BFGS:   2  16:14:26        0.620874      13.0140
  BFGS:   3  16:14:26       -0.028619       4.4019
  BFGS:   4  16:14:26       -0.129349       0.7307
  BFGS:   5  16:14:26       -0.132320       0.0138

When doing structure optimization, it is useful to write the
trajectory to a file, so that the progress of the optimization run can
be followed during or after the run::

  dyn = BFGS(water, trajectory='H2O.traj')
  dyn.run(fmax=0.05)
  
Use the command ``ag H2O.traj`` to see what is going on (more here:
:mod:`gui`).  The trajectory file can also be accessed using the
module :mod:`ase.io.trajectory`.

The ``attach`` method takes an optional argument ``interval=n`` that can
be used to tell the structure optimizer object to write the
configuration to the trajectory file only every ``n`` steps.

During a structure optimization, the :class:`BFGS` and
:class:`LBFGS` optimizers use two quantities to decide where to move
the atoms on each step:

 * the forces on each atom, as returned by the associated :class:`Calculator`
   object
 * the Hessian matrix, i.e. the matrix of second derivatives
   :math:`\frac{\partial^2 E}{\partial x_i \partial x_j}` of the
   total energy with respect to nuclear coordinates.

If the atoms are close to the minimum, such that the potential energy
surface is locally quadratic, the Hessian and forces accurately
determine the required step to reach the optimal structure.  The
Hessian is very expensive to calculate *a priori*, so instead the
algorithm estimates it by means of an initial guess which is adjusted
along the way depending on the information obtained on each step of
the structure optimization.

It is frequently practical to restart or continue a structure
optimization with a geometry obtained from a previous relaxation.
Aside from the geometry, the Hessian of the previous run can and
should be retained for the second run.  Use the ``restart`` keyword to
specify a file in which to save the Hessian::

  dyn = BFGS(system, trajectory='qn.traj', restart='qn.pckl')

This will create an optimizer which saves the Hessian to
:file:`qn.pckl` (using the Python :mod:`pickle` module) on each
step.  If the file already exists, the Hessian will also be
*initialized* from that file.

The trajectory file can also be used to restart a structure
optimization, since it contains the history of all forces and
positions, and thus whichever information about the Hessian was
assembled so far::

  dyn = BFGS(system, trajectory='qn.traj')
  dyn.replay_trajectory('history.traj')

This will read through each iteration stored in :file:`history.traj`,
performing adjustments to the Hessian as appropriate.  Note that these
steps will not be written to :file:`qn.traj`.  If restarting with more than
one previous trajectory file, use :command:`ag` to concatenate them
into a single trajectory file first::

  $ ag part1.traj part2.traj -o history.traj

The file :file:`history.traj` will then contain all necessary
information.

When switching between different types of optimizers, e.g. between
``BFGS`` and ``LBFGS``, the pickle-files specified by the
``restart`` keyword are not compatible, but the Hessian can still be
retained by replaying the trajectory as above.

.. note::

   In many of the examples, tests, exercises and tutorials,
   ``QuasiNewton`` is used -- it is a synonym for ``BFGS``.


LBFGS
-----
.. module:: optimize.lbfgs

LBFGS is the limited memory version of the BFGS algorithm, where 
the inverse of Hessian matrix is updated instead of the Hessian
itself. Two ways exist for determining the atomic
step: Standard ``LBFGS`` and ``LBFGSLineSearch``. For the 
first one, both the directions and lengths of the atomic steps 
are determined by the approximated Hessian matrix. While for the 
latter one, the approximated Hessian matrix is only used to find 
out the directions of the line searches and atomic steps, the 
step lengths are determined by the forces. 

To start a structure optimization with LBFGS algorithm is similar to
BFGS. A typical optimization should look like::

  dyn = LBFGS(system, trajectory='lbfgs.traj', restart='lbfgs.pckl')

where the trajectory and the restart save the trajectory of the 
optimization and the vectors needed to generate the Hessian Matrix.


FIRE
----
.. module:: optimize.fire

Read about this algorithm here:

  | Erik Bitzek, Pekka Koskinen, Franz Gähler, Michael Moseler, and Peter Gumbsch
  | `Structural Relaxation Made Simple`__
  | Physical Review Letters, Vol. **97**, 170201 (2006)

__ http://dx.doi.org/10.1103/PhysRevLett.97.170201


MDMin
-----
.. module:: optimize.mdmin

The MDmin algorithm is a modification of the usual velocity-Verlet
molecular dynamics algorithm.  Newtons second law is solved
numerically, but after each time step the dot product between the
forces and the momenta is checked.  If it is zero, the system has just
passed through a (local) minimum in the potential energy, the kinetic
energy is large and about to decrease again.  At this point, the
momentum is set to zero.  Unlike a "real" molecular dynamics, the
masses of the atoms are not used, instead all masses are set to one.

The MDmin algorithm exists in two flavors, one where each atom is
tested and stopped individually, and one where all coordinates are
treated as one long vector, and all momenta are set to zero if the
dotproduct between the momentum vector and force vector (both of
length 3N) is zero.  This module implements the latter version.

Although the algorithm is primitive, it performs very well because it
takes advantage of the physics of the problem.  Once the system is so
near the minimum that the potential energy surface is approximately
quadratic it becomes advantageous to switch to a minimization method
with quadratic convergence, such as `Conjugate Gradient` or `Quasi
Newton`.


SciPy optimizers
----------------
.. module:: optimize.sciopt

SciPy provides a number of optimizers. An interface module for a couple of
these have been written for ASE. Most notable are the optimizers SciPyFminBFGS
and SciPyFminCG. These are called with the regular syntax and can be imported
as::

  from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG

.. autoclass:: ase.optimize.sciopt.SciPyFminBFGS
.. autoclass:: ase.optimize.sciopt.SciPyFminCG

.. seealso::

  :epydoc:`optimize.sciopt.SciPyFminBFGS`, 
  :epydoc:`optimize.sciopt.SciPyFminCG`


BFGSLineSearch
--------------
.. module:: optimize.bfgslinesearch

BFGSLineSearch is the BFGS algorithm with an line search mechanism
that enforces the step taken fulfills the Wolfe conditions, so that
the energy and absolute value of the force decrease monotonically. Like
the lbfgs algorithm the inverse of the Hessian Matrix is updated.

The usage of BFGSLineSearch algorithm is similar to other BFGS type
algorithms. A typical optimization should look like::

  from ase.optimize.bfgslinesearch import BFGSLineSearch

  dyn = BFGSLineSearch(system, trajectory='bfgs_ls.traj', restart='bfgs_ls.pckl')

where the trajectory and the restart save the trajectory of the
optimization and the information needed to generate the Hessian Matrix.


Global optimization
===================

There is currently one global optimisation algorithm available.


Basin hopping
-------------
.. module:: optimize.basin

The global optimization algorithm can be used quite similar as a 
local optimization algorithm::

  from ase import *
  from ase.optimize.basin import BasinHopping

  bh = BasinHopping(system,               # the system to optimize 
                    temperature=100 * kB, # 'temperature' to overcome barriers
                    dr=0.5,               # maximal stepwidth
	       	    optimizer=LBFGS,      # optimizer to find local minima
		    fmax=0.1,             # maximal force for the optimizer
                    )

Read more about this algorithm here:

  | David J. Wales and Jonathan P. K. Doye
  | `Global Optimization by Basin-Hopping and the Lowest Energy Structures of Lennard-Jones Clusters Containing up to 110 Atoms`__
  | J. Phys. Chem. A, Vol. **101**, 5111-5116 (1997)

__ http://pubs.acs.org/doi/abs/10.1021/jp970984n

and here:

  | David J. Wales and Harold A. Scheraga
  | `Global Optimization of Clusters, Crystals, and Biomolecules`__
  | Science, Vol. **285**, 1368 (1999)

__ http://www.sciencemag.org/cgi/content/abstract/sci;285/5432/1368


