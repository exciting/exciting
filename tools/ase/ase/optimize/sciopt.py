import numpy as np
try:
    import scipy.optimize as opt
except ImportError:
    pass

from ase.optimize.optimize import Optimizer


class Converged(Exception):
    pass

class OptimizerConvergenceError(Exception):
    pass

class SciPyOptimizer(Optimizer):
    """General interface for SciPy optimizers

    Only the call to the optimizer is still needed
    """
    def __init__(self, atoms, logfile='-', trajectory=None,
                 callback_always=False, alpha=70.0):
        """Initialize object

        Parameters:

        callback_always: book
            Should the callback be run after each force call (also in the
            linesearch)

        alpha: float
            Initial guess for the Hessian (curvature of energy surface). A
            conservative value of 70.0 is the default, but number of needed
            steps to converge might be less if a lower value is used. However,
            a lower value also means risk of instability.

        """
        restart = None
        Optimizer.__init__(self, atoms, restart, logfile, trajectory)
        self.force_calls = 0
        self.callback_always = callback_always
        self.H0 = alpha

    def x0(self):
        """Return x0 in a way SciPy can use

        This class is mostly usable for subclasses wanting to redefine the
        parameters (and the objective function)"""
        return self.atoms.get_positions().reshape(-1)

    def f(self, x):
        """Objective function for use of the optimizers"""
        self.atoms.set_positions(x.reshape(-1, 3))
        # Scale the problem as SciPy uses I as initial Hessian.
        return self.atoms.get_potential_energy() / self.H0

    def fprime(self, x):
        """Gradient of the objective function for use of the optimizers"""
        self.atoms.set_positions(x.reshape(-1, 3))
        self.force_calls += 1

        if self.callback_always:
            self.callback(x)

        # Remember that forces are minus the gradient!
        # Scale the problem as SciPy uses I as initial Hessian.
        return - self.atoms.get_forces().reshape(-1) / self.H0

    def callback(self, x):
        """Callback function to be run after each iteration by SciPy

        This should also be called once before optimization starts, as SciPy
        optimizers only calls it after each iteration, while ase optimizers
        call something similar before as well.
        """
        f = self.atoms.get_forces()
        self.log(f)
        self.call_observers()
        if self.converged(f):
            raise Converged
        self.nsteps += 1

    def run(self, fmax=0.05, steps=100000000):
        self.fmax = fmax
        # As SciPy does not log the zeroth iteration, we do that manually
        self.callback(None)
        try:
            # Scale the problem as SciPy uses I as initial Hessian.
            self.call_fmin(fmax / self.H0, steps)
        except Converged:
            pass

    def dump(self, data):
        pass

    def load(self):
        pass

    def call_fmin(self, fmax, steps):
        raise NotImplementedError

class SciPyFminCG(SciPyOptimizer):
    """Non-linear (Polak-Ribiere) conjugate gradient algorithm"""
    def call_fmin(self, fmax, steps):
        output = opt.fmin_cg(self.f,
                             self.x0(),
                             fprime=self.fprime,
                             #args=(),
                             gtol=fmax * 0.1, #Should never be reached
                             norm=np.inf,
                             #epsilon=
                             maxiter=steps,
                             full_output=1,
                             disp=0,
                             #retall=0, 
                             callback=self.callback
                            )
        warnflag = output[-1]
        if warnflag == 2:
            raise OptimizerConvergenceError('Warning: Desired error not necessarily achieved ' \
                                            'due to precision loss')

class SciPyFminBFGS(SciPyOptimizer):
    """Quasi-Newton method (Broydon-Fletcher-Goldfarb-Shanno)"""
    def call_fmin(self, fmax, steps):
        output = opt.fmin_bfgs(self.f,
                               self.x0(),
                               fprime=self.fprime,
                               #args=(), 
                               gtol=fmax * 0.1, #Should never be reached
                               norm=np.inf,
                               #epsilon=1.4901161193847656e-08, 
                               maxiter=steps,
                               full_output=1,
                               disp=0,
                               #retall=0, 
                               callback=self.callback
                              )
        warnflag = output[-1]
        if warnflag == 2:
            raise OptimizerConvergenceError('Warning: Desired error not necessarily achieved' \
                                            'due to precision loss')

class SciPyGradientlessOptimizer(Optimizer):
    """General interface for gradient less SciPy optimizers

    Only the call to the optimizer is still needed

    Note: If you redefien x0() and f(), you don't even need an atoms object.
    Redefining these also allows you to specify an arbitrary objective
    function.

    XXX: This is still a work in progress
    """
    def __init__(self, atoms, logfile='-', trajectory=None,
                 callback_always=False):
        """Parameters:

        callback_always: book
            Should the callback be run after each force call (also in the
            linesearch)
        """
        restart = None
        Optimizer.__init__(self, atoms, restart, logfile, trajectory)
        self.function_calls = 0
        self.callback_always = callback_always

    def x0(self):
        """Return x0 in a way SciPy can use

        This class is mostly usable for subclasses wanting to redefine the
        parameters (and the objective function)"""
        return self.atoms.get_positions().reshape(-1)

    def f(self, x):
        """Objective function for use of the optimizers"""
        self.atoms.set_positions(x.reshape(-1, 3))
        self.function_calls += 1
        # Scale the problem as SciPy uses I as initial Hessian.
        return self.atoms.get_potential_energy()

    def callback(self, x):
        """Callback function to be run after each iteration by SciPy

        This should also be called once before optimization starts, as SciPy
        optimizers only calls it after each iteration, while ase optimizers
        call something similar before as well.
        """
        # We can't assume that forces are available!
        #f = self.atoms.get_forces()
        #self.log(f)
        self.call_observers()
        #if self.converged(f):
        #    raise Converged
        self.nsteps += 1

    def run(self, ftol=0.01, xtol=0.01, steps=100000000):
        self.xtol = xtol
        self.ftol = ftol
        # As SciPy does not log the zeroth iteration, we do that manually
        self.callback(None)
        try:
            # Scale the problem as SciPy uses I as initial Hessian.
            self.call_fmin(xtol, ftol, steps)
        except Converged:
            pass

    def dump(self, data):
        pass

    def load(self):
        pass

    def call_fmin(self, fmax, steps):
        raise NotImplementedError

class SciPyFmin(SciPyGradientlessOptimizer):
    """Nelder-Mead Simplex algorithm

    Uses only function calls.

    XXX: This is still a work in progress
    """
    def call_fmin(self, xtol, ftol, steps):
        output = opt.fmin(self.f,
                          self.x0(),
                          #args=(),
                          xtol=xtol,
                          ftol=ftol,
                          maxiter=steps,
                          #maxfun=None,
                          #full_output=1,
                          disp=0,
                          #retall=0,
                          callback=self.callback
                         )

class SciPyFminPowell(SciPyGradientlessOptimizer):
    """Powell's (modified) level set method

    Uses only function calls.

    XXX: This is still a work in progress
    """
    def __init__(self, *args, **kwargs):
        """Parameters:

        direc: float
            How much to change x to initially. Defaults to 0.04.
        """
        direc = kwargs.pop('direc', None)
        SciPyGradientlessOptimizer.__init__(self, *args, **kwargs)

        if direc is None:
            self.direc = np.eye(len(self.x0()), dtype=float) * 0.04
        else:
            self.direc = np.eye(len(self.x0()), dtype=float) * direc

    def call_fmin(self, xtol, ftol, steps):
        output = opt.fmin_powell(self.f,
                                 self.x0(),
                                 #args=(),
                                 xtol=xtol,
                                 ftol=ftol,
                                 maxiter=steps,
                                 #maxfun=None,
                                 #full_output=1,
                                 disp=0,
                                 #retall=0,
                                 callback=self.callback,
                                 direc=self.direc
                                )
