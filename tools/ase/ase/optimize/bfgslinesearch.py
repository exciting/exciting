#__docformat__ = "restructuredtext en"
# ******NOTICE***************
# optimize.py module by Travis E. Oliphant
#
# You may copy and use this module as you see fit with no
# guarantee implied provided you keep this notice in all copies.
# *****END NOTICE************

import numpy as np
from numpy import atleast_1d, eye, mgrid, argmin, zeros, shape, empty, \
     squeeze, vectorize, asarray, absolute, sqrt, Inf, asfarray, isinf
from ase.utils.linesearch import LineSearch
from ase.optimize.optimize import Optimizer
from numpy import arange


# These have been copied from Numeric's MLab.py
# I don't think they made the transition to scipy_core

# Modified from scipy_optimize
abs = absolute
import __builtin__
pymin = __builtin__.min
pymax = __builtin__.max
__version__="0.1"

class BFGSLineSearch(Optimizer):
    def __init__(self, atoms, restart=None, logfile='-', maxstep=.2,
                 trajectory=None, c1=.23, c2=0.46, alpha=10., stpmax=50.,
                 use_free_energy=True):
        """Minimize a function using the BFGS algorithm.

        Notes:

            Optimize the function, f, whose gradient is given by fprime
            using the quasi-Newton method of Broyden, Fletcher, Goldfarb,
            and Shanno (BFGS) See Wright, and Nocedal 'Numerical
            Optimization', 1999, pg. 198.

        *See Also*:

          scikits.openopt : SciKit which offers a unified syntax to call
                            this and other solvers.

        """
        self.maxstep = maxstep
        self.stpmax = stpmax
        self.alpha = alpha
        self.H = None
        self.c1 = c1
        self.c2 = c2
        self.force_calls = 0
        self.function_calls = 0
        self.r0 = None
        self.g0 = None
        self.e0 = None
        self.load_restart = False
        self.task = 'START'
        self.rep_count = 0
        self.p = None
        self.alpha_k = None
        self.no_update = False
        self.replay = False
        self.use_free_energy = use_free_energy

        Optimizer.__init__(self, atoms, restart, logfile, trajectory)

    def read(self):
        self.r0, self.g0, self.e0, self.task, self.H = self.load()
        self.load_restart = True    

    def reset(self):
        print 'reset'
        self.H = None
        self.r0 = None
        self.g0 = None
        self.e0 = None
        self.rep_count = 0
          

    def step(self, f):
        atoms = self.atoms
        from ase.neb import NEB 
        assert not isinstance(atoms, NEB) 
        r = atoms.get_positions()
        r = r.reshape(-1)
        g = -f.reshape(-1) / self.alpha
        p0 = self.p
        self.update(r, g, self.r0, self.g0, p0)
        #o,v = np.linalg.eigh(self.B)
        e = self.func(r)

        self.p = -np.dot(self.H,g)
        p_size = np.sqrt((self.p **2).sum())
        if self.nsteps != 0:
            p0_size = np.sqrt((p0 **2).sum())
            delta_p = self.p/p_size + p0/p0_size
        if p_size <= np.sqrt(len(atoms) * 1e-10):
            self.p /= (p_size / np.sqrt(len(atoms)*1e-10))
        ls = LineSearch()
        self.alpha_k, e, self.e0, self.no_update = \
           ls._line_search(self.func, self.fprime, r, self.p, g, e, self.e0,
                           maxstep=self.maxstep, c1=self.c1,
                           c2=self.c2, stpmax=self.stpmax)
        if self.alpha_k is None:
            raise RuntimeError("LineSearch failed!")

        dr = self.alpha_k * self.p
        atoms.set_positions((r+dr).reshape(len(atoms),-1))
        self.r0 = r
        self.g0 = g
        self.dump((self.r0, self.g0, self.e0, self.task, self.H))

    def update(self, r, g, r0, g0, p0):
        self.I = eye(len(self.atoms) * 3, dtype=int)
        if self.H is None:
            self.H = eye(3 * len(self.atoms))
            #self.B = np.linalg.inv(self.H)
            return
        else:
            dr = r - r0
            dg = g - g0 
            if not ((self.alpha_k > 0 and abs(np.dot(g,p0))-abs(np.dot(g0,p0)) < 0) \
                or self.replay):
                return
            if self.no_update == True:
                print 'skip update'
                return

            try: # this was handled in numeric, let it remaines for more safety
                rhok = 1.0 / (np.dot(dg,dr))
            except ZeroDivisionError:
                rhok = 1000.0
                print "Divide-by-zero encountered: rhok assumed large"
            if isinf(rhok): # this is patch for np
                rhok = 1000.0
                print "Divide-by-zero encountered: rhok assumed large"
            A1 = self.I - dr[:, np.newaxis] * dg[np.newaxis, :] * rhok
            A2 = self.I - dg[:, np.newaxis] * dr[np.newaxis, :] * rhok
            H0 = self.H
            self.H = np.dot(A1, np.dot(self.H, A2)) + rhok * dr[:, np.newaxis] \
                     * dr[np.newaxis, :]
            #self.B = np.linalg.inv(self.H)

    def func(self, x):
        """Objective function for use of the optimizers"""
        self.atoms.set_positions(x.reshape(-1, 3))
        calc = self.atoms.get_calculator()
        self.function_calls += 1
        # Scale the problem as SciPy uses I as initial Hessian.
        if self.use_free_energy:
            try:
                return calc.get_potential_energy(self.atoms,force_consistent=True) / self.alpha
            except TypeError:
                return calc.get_potential_energy(self.atoms) / self.alpha
        else:
            return calc.get_potential_energy(self.atoms) / self.alpha
    
    def fprime(self, x):
        """Gradient of the objective function for use of the optimizers"""
        self.atoms.set_positions(x.reshape(-1, 3))
        self.force_calls += 1
        # Remember that forces are minus the gradient!
        # Scale the problem as SciPy uses I as initial Hessian.
        f = self.atoms.get_forces().reshape(-1) 
        return - f / self.alpha

    def replay_trajectory(self, traj):
        """Initialize hessian from old trajectory."""
        self.replay = True
        if isinstance(traj, str):
            from ase.io.trajectory import PickleTrajectory
            traj = PickleTrajectory(traj, 'r')
        atoms = traj[0]
        r0 = None
        g0 = None
        for i in range(0, len(traj) - 1):
            r = traj[i].get_positions().ravel()
            g = - traj[i].get_forces().ravel() / self.alpha
            self.update(r, g, r0, g0, self.p)
            self.p = -np.dot(self.H,g)
            r0 = r.copy()
            g0 = g.copy()
        self.r0 = r0
        self.g0 = g0

def wrap_function(function, args):
    ncalls = [0]
    def function_wrapper(x):
        ncalls[0] += 1
        return function(x, *args)
    return ncalls, function_wrapper
