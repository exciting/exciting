.. _lattice_constant:

=========================
Finding lattice constants
=========================

.. seealso::

   :ref:`eos`.


HCP
===

Let's try to find the `a` and `c` lattice constants for HCP nickel
using the :mod:`EMT <emt>` potential.  

First, we make a good intial guess for `a` and `c` using the FCC nearest
neighbor distance and the ideal `c/a` ratio::

  from numpy import sqrt
  a0 = 3.52 / sqrt(2)
  c0 = sqrt(8 / 3.0) * a0

and create a trajectory for the results::

  from ase.io import PickleTrajectory
  traj = PickleTrajectory('Ni.traj', 'w')

Finally, we do the 12 calculations (four values for `a` and three for `c`)::

  import numpy as np
  from ase.structure import bulk
  from ase.calculators import EMT
  eps = 0.01
  for a in a0 * np.linspace(1 - eps, 1 + eps, 4):
      for c in c0 * np.linspace(1 - eps, 1 + eps, 3):
          ni = bulk('Ni', 'hcp', a=a, covera=c / a)
          ni.set_calculator(EMT())
          ni.get_potential_energy()
          traj.write(ni)


Analysis
--------

Now, we need to extract the data from the trajectory.  Try this:

>>> from ase.structure import bulk
>>> ni = bulk('Ni', 'hcp', a=2.5, covera=4.0 / 2.5)
>>> ni.cell
array([[ 2.5       ,  0.        ,  0.        ],
       [ 1.25      ,  2.16506351,  0.        ],
       [ 0.        ,  0.        ,  4.        ]])

So, we can get `a` and `c` from ``ni.cell[0, 0]`` and ``ni.cell[2,
2]``:

>>> from ase.io import read
>>> configs = read('Ni.traj@:')
>>> energies = [config.get_potential_energy() for config in configs]
>>> ac = [(config.cell[0, 0], config.cell[2, 2]) for config in configs]

We fit the energy to this expression:

.. math:: c_0 + c_1 a + c_2 c + c_3 a^2 + c_4 ac + c_5 c^2 +
          c_6 a^3 + c_7 a^2c + c_8 ac^2 + c_9 c^3

>>> from ase.optimize import polyfit
>>> p = polyfit(ac, energies)

using the function:

.. autofunction:: ase.optimize.polyfit

The minimum can be found using SciPy's fmin_bfgs_
function:

>>> from scipy.optimize import fmin_bfgs
>>> a0 = 3.52 / sqrt(2)
>>> c0 = sqrt(8 / 3.0) * a0
>>> a, c = fmin_bfgs(p, (a0, c0))
Warning: Desired error not necessarily achieveddue to precision loss
         Current function value: 0.010030
         Iterations: 7
         Function evaluations: 425
         Gradient evaluations: 106
>>> print a, c
2.46888950503 4.02027198125

In (often) cases optimization fails, one may use
another energy expression:

.. math:: c_0 + c_1 a + c_2 c + c_3 a^2 + c_4 c^2

>>> import numpy as np
>>> sle = np.linalg.solve
>>> E = np.array(energies)
>>> A = np.array([(1, x, y, x**2, y**2)
>>>               for x, y in ac]).T
>>> C = sle(np.inner(A, A), np.dot(A, E))

>>> a = - C[1] / (2 * C[3])
>>> c = - C[2] / (2 * C[4])
>>> print a, c
2.46465936501 4.00438337976

.. _fmin_bfgs: http://docs.scipy.org/doc/scipy/reference/optimize.html
