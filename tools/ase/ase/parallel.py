import sys
import time
import atexit

import numpy as np


def paropen(name, mode='r', buffering=0):
    """MPI-safe version of open function.

    In read mode, the file is opened on all nodes.  In write and
    append mode, the file is opened on the master only, and /dev/null
    is opened on all other nodes.
    """
    if rank > 0 and mode[0] != 'r':
        name = '/dev/null'
    return open(name, mode, buffering)


def parprint(*args, **kwargs):
    """MPI-safe print - prints only from master.

    Tries to adopt python 3 behaviour.
    """
    if rank > 0:
        return
    defaults = {'end': '\n',
                'file': sys.stdout }
    for key in defaults:
        if not key in kwargs:
            kwargs[key] = defaults[key]

    for arg in args[:-1]:
        print >> kwargs['file'], arg,
    if len(args):
        last = args[-1]
    else:
        last = ''
    if kwargs['end'] == '\n':
        print >> kwargs['file'], last
    else:
        print >> kwargs['file'], last,


class DummyMPI:
    rank = 0
    size = 1
    def sum(self, a):
        if isinstance(a, np.ndarray) and a.ndim > 0:
            pass
        else:
            return a
    
    def barrier(self):
        pass

    def broadcast(self, a, rank):
        pass


class MPI4PY:
    def __init__(self):
        from mpi4py import MPI
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.rank
        self.size = self.comm.size

    def sum(self, a):
        return self.comm.allreduce(a)
    
    def barrier(self):
        self.comm.barrier()

    def abort(self, code):
        self.comm.Abort(code)

    def broadcast(self, a, rank):
        a[:] = self.comm.bcast(a, rank)


# Check for special MPI-enabled Python interpreters:
if '_gpaw' in sys.modules:
    # http://wiki.fysik.dtu.dk/gpaw
    from gpaw.mpi import world
elif 'asapparallel3' in sys.modules:
    # http://wiki.fysik.dtu.dk/Asap
    # We cannot import asap3.mpi here, as that creates an import deadlock
    #from asap3.mpi import world
    import asapparallel3
    world = asapparallel3.Communicator()
elif 'Scientific_mpi' in sys.modules:
    from Scientific.MPI import world
elif 'mpi4py' in sys.modules:
    world = MPI4PY()
else:
    # This is a standard Python interpreter:
    world = DummyMPI()

rank = world.rank
size = world.size
barrier = world.barrier


def register_parallel_cleanup_function():
    """Call MPI_Abort if python crashes.

    This will terminate the processes on the other nodes."""
        
    if size == 1:
        return

    def cleanup(sys=sys, time=time, world=world):
        error = getattr(sys, 'last_type', None)
        if error:
            sys.stdout.flush()
            sys.stderr.write(('ASE CLEANUP (node %d): %s occurred.  ' +
                              'Calling MPI_Abort!\n') % (world.rank, error))
            sys.stderr.flush()
            # Give other nodes a moment to crash by themselves (perhaps
            # producing helpful error messages):
            time.sleep(3)
            world.abort(42)

    atexit.register(cleanup)

def distribute_cpus(parsize_calculator, comm):
    """Distribute cpus to tasks and calculators.

    parsize_calculator: number of nodes per calculator

    comm: total communicator object
    """
    
    assert parsize_calculator <= comm.size
    assert comm.size % parsize_calculator == 0

    tasks_rank = comm.rank // parsize_calculator

    r0 = tasks_rank * parsize_calculator
    ranks = np.arange(r0, r0 + parsize_calculator)
    calc_comm = comm.new_communicator(ranks)

    tasks_comm = np.arange(0, comm.size, parsize_calculator)

    return calc_comm, tasks_comm, tasks_rank
