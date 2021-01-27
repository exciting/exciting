"""
Parser command-line arguments and run the test suite.

To run all tests in serial, type `python3 runtest.py -a run`
To run all tests in parallel, type `python3 runtest.py -a run -e excitingmpismp`
For more details, type `python3 runtest.py --help`
"""

import sys
import os
import argparse as ap
import unittest
import warnings
from collections import namedtuple

sys.path.insert(1, 'tools')
sys.path.insert(1, 'tools/selftests/')

from termcolor_wrapper import print_color
from procedures import *
from constants import settings
import runselftests
from failing_tests import skipped_tests  


def optionParser(testFarm:str, exedir:str):
    """
    Parse command line inputs 
    
    Parse: 
      action to perform
      test/s to run
      executable choice
      number of MPI processes
      number of OMP threads

    :param testFarm: test farm directory
    :param exedir: exciting executable directory 
    :return input_options: Dictionary of parsed command line arguments 
    """
 
    p = ap.ArgumentParser(description="Usage: python3 runtest.py -a <action> -t <tests> -e <execuatble> -np <NP> -omp <omp>")

    help_action = "Defines what action is done. " \
        + "'run' for running tests; " \
        + "'ref' for running references; " \
        + "'clean' for cleaning the test directories; " \
        + "'report' for collecting results; default is 'run'"

    p.add_argument ('-a',
                    metavar = '--action',
                    help = help_action,
                    type = str,
                    choices = settings.action_choices,
                    default = 'run')

    help_test = "Test cases for the action. This can be either: " \
        + "<test case 1>, <test case 2>, ... for special test cases. " \
        + "This are the names of the test case directories without '/'; " \
        + "'all' for all test cases; default is 'all'"
    
    p.add_argument('-t',
                   metavar ='--tests',
                   help = help_test,
                   type = str,
                   default = ['all'] ,
                   nargs = '+')

    help_executable = "exciting executables. "\
        + "'excitingser' for serial version of the code; "\
        + "'excitingmpi' for version with MPI parallisation; "\
        + "'excitingmpismp for version with MPI amd SMP parallisation;" \
        + "Default is excitingser"
    
    p.add_argument('-e',
                   metavar = '--executable',
                   help = help_executable,
                   type = str,
                   default = 'excitingser',
                   choices = settings.binary_names)
    
    p.add_argument('-np',
                   metavar = '--NP',
                   help = "Number of cores for MPI run. Can only be used in " + \
                   "combination with excitingmpi or excitingmpismp as executable. " + \
                   "Default is 2 for MPI and MPI+OMP calculations, and 1 for serial or pure OMP", 
                   type = int)
    
    p.add_argument('-omp', metavar='--ompthreads', help='Number of OMP threads. ' + \
                   'Default is 2 for OMP and MPI+OMP calculations, and 1 for serial or pure MPI',
                   type=int)

    args = p.parse_args()

    input_options = {'action': args.a, 'tests':args.t}
    input_options['np']  = args.np if args.np is not None else settings.default_np[args.e] 
    input_options['omp'] = args.omp if args.omp is not None else settings.default_threads[args.e]

    if args.e == 'excitingser':
        if args.np is not None:
            warnings.warn('Using serial exciting, -np will be ignored.')
        if args.omp is not None:
            warnings.warn('Using serial exciting, -omp will be ignored.')
    if args.e == 'excitingmpi':
         if args.omp is not None:
            warnings.warn('Using pure mpi exciting, -omp will be ignored.')
            
    if input_options['tests'][0] == 'all':
        input_options['tests'] = next(os.walk(testFarm))[1]
    else:
        input_options['tests'] = input_options['tests']
        
    testDirs = next(os.walk(testFarm))[1]
    for t in input_options['tests']:
        if t not in testDirs:
            raise ValueError('Could not find %s/ in %s/'%(t, testFarm))

    input_options['executable'] = os.path.join(exedir, args.e)            
    if args.e in ['excitingmpi', 'excitingmpismp']:
        input_options['executable'] = \
            'mpirun -np %i %s' % (input_options['np'], input_options['executable'])

    return input_options 


def main(settings:namedtuple, input_options:dict):
    """
    Run test suite actions
 
    These include 'run' the test suite, rerun 'ref'erence calculations
    or collect 'report's, depending on the action specified as a 
    command-line argument.   

    :param settings: Default settings for environment variables, paths
                     and run-time parameters. 
    :param input_options: parsed command-line arguments
    """

    print('Run test suite:')

    os.environ["OMP_NUM_THREADS"] = str(input_options['omp'])
    
    species_files = next(os.walk(settings.species))[2]
    action = input_options['action']
    executable_command = input_options['executable']
    executable = executable_command.split('/')[-1]

    if action == "run":
        runTests(settings.test_farm, 
                 settings.main_output,
                 input_options['tests'],
                 settings.run_dir, 
                 settings.ref_dir,
                 settings.init_default,
                 executable_command,
                 input_options['np'],
                 input_options['omp'],
                 settings.max_time,
                 skipped_tests[executable])
     
    elif action == "ref":    
        while(True):
            answer = input("Are you sure that you want to rerun reference(s)? " + \
                           "This will replace the old reference(s) with reference(s) " + \
                           "created by the current state of the code and should only " + \
                           "be done if all tests pass with the current version of the code.(y/n)")
            if answer=='y':
                break
            elif answer=='n':
                sys.exit()
        runReferences(settings.test_farm,
                      settings.main_output,
                      input_options['tests'],
                      settings.ref_dir,
                      executable_command,
                      settings.ignored_output + species_files,
                      input_options['np'],
                      input_options['omp'],
                      settings.max_time)
        collectReports(settings.test_farm, input_options['tests'])
       
    elif action == "report":
        collectReports(settings.test_farm, input_options['tests'])
    

def run_self_tests(verbosity=3):
    """
    Test the application test suite functions 
    Execute selftests before the test suite
   
    :return: result  
    """
    print('Run self tests:')
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add here more self tests
    suite.addTests(loader.loadTestsFromModule(runselftests))
    
    runner = unittest.TextTestRunner(verbosity=verbosity)
    result = runner.run(suite)
    
    if len(result.errors)==0:
        print_color('Self tests SUCCESS', 'green')
    else:
        print_color('Self tests FAIL', 'red')
        quit('Test suite has quit')
        
    return result


if __name__ == "__main__":

    input_options = optionParser(settings.test_farm, settings.exe_dir)

    # make clean
    if input_options['action'] == "clean":
        species_files = next(os.walk(settings.species))[2]
        
        cleanTests(settings.test_farm, 
                   input_options['tests'],
                   settings.run_dir, 
                   settings.ref_dir,
                   settings.not_clean + species_files)
        exit()
        
    self_test_results = run_self_tests()
    
    main(settings, input_options)
