""" Execute a program.
"""
import os
from subprocess import PIPE, Popen, TimeoutExpired
import time
from typing import Tuple

from ..exciting_settings.constants import RunProperties


def set_job_environment(threads: dict):
    """ Create an environment instance.

    :param threads: OMP and MKL threads.
    :return Instance of os._Environ (behaves like a dict)
    """
    my_env = os.environ.copy()

    my_env["OMP_NUM_THREADS"] = str(threads['omp'])

    mkl_threads = threads['mkl_threads']
    if mkl_threads is not None:
        my_env["MKL_NUM_THREADS"] = str(mkl_threads)

    return my_env


def execute_job(job: RunProperties, my_env=None) -> Tuple[bool, str, float]:
    """Executes a calculation and checks if it was successful.

    :param job: Run properties.
    :param my_env: Optional environment instance.
    :return (run_success, err_mess, run_time): Successful calculation, stdout and run time.
    """
    # Use existing environment
    if my_env is None:
        my_env = os.environ.copy()

    terminated_cleanly, err_mess, run_time = execute(job.run_dir, job.executable_cmd, job.max_time, my_env)
    run_success = job.calculation_completed(terminated_cleanly)

    return run_success, err_mess, run_time


def execute(path, execution_str: str, max_time: int, my_env) -> Tuple[bool, str, float]:
    """Executes a calculation run and checks if it terminated cleanly.

    :param path: Run directory.
    :param execution_str: Execution string. For example:
    `binary.exe` or `mpirun -np NP binary.exe`.
    :param max_time: Maximum time in seconds after which the calculation is considered to have timed out.
    :param my_env: Shell environment instance.
    :return (terminated_cleanly, err_mess, run_time): Successful calculation, stdout and run time.
    """
    t_start = time.time()

    process = Popen(execution_str.split(), cwd=path, stdout=PIPE, stderr=PIPE, env=my_env)

    try:
        err_mess = process.communicate(timeout=max_time)[0]
    except TimeoutExpired:
        process.kill()
        err_mess = 'Time expired.'

    process.wait()
    terminated_cleanly = process.returncode == 0
    t_end = time.time()

    return terminated_cleanly, err_mess, t_end - t_start
