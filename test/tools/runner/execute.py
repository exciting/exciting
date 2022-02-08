#TODO(Bene): Robust method to check if a calculation was successfull or not. 
#                So far the functions checks if INFO.OUT is present but that is not good
#                since the run can crash in a scl cycle.

import os
from subprocess import PIPE, Popen, TimeoutExpired
import time


def execute(path:str, executable_str:str, mainOut:str, maxTime:int):
    """
    Executes a exciting run, checks if it was successfull.
    :param path:          path to directory where the executable will be run
    :param executable_str:    executable command (poorly-named). For example:
                              exciting_serial exciting_smp or mpirun -np NP exciting_mpismp
    :param mainOut:       main output file (INFO.OUT)
    :param maxTime:       maximum time for an exciting run in seconds
    
    Output:
        success     bool                true if run was successfull, false else
        errMess     list of strings     terminal output of exciting
        run_time    float     Run time of job 
    """
    current_dir = os.getcwd()

    try:
        os.chdir(os.path.join(current_dir, path))
    except OSError:
        raise OSError('Could not enter %s.'%path)

    t_start = time.time()
    executable = Popen(executable_str.split(), stdout=PIPE)
    
    try:
        errMess = executable.communicate(timeout=maxTime)[0]
        runSucc = os.path.isfile(mainOut)
    except TimeoutExpired:
        executable.kill()
        errMess = 'Time expired.'
        runSucc = False
    
    executable.wait()
    t_end = time.time()
    os.chdir(current_dir)
    
    return runSucc, errMess, t_end-t_start