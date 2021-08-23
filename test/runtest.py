"""
Parser command-line arguments and run the test suite.

To run all tests in serial, type `python3 runtest.py -a run`
To run all tests in parallel, type `python3 runtest.py -a run -e exciting_mpismp`
For more details, type `python3 runtest.py --help`
"""
import sys
import argparse as ap
import unittest
import warnings
from collections import namedtuple
import os


from tools.termcolor_wrapper import print_color
from tools.runner.test import run_tests
from tools.runner.reference import run_references
from tools.runner.clean import clean_tests
from tools.constants import settings
from tools.selftests import runselftests
from tools.parsers import install_excitingtools
from tools.utils import Build_type, build_type_str_to_enum, build_type_enum_to_str
from failing_tests import set_skipped_tests

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return '%s: %s \n' % (category.__name__, message)

warnings.formatwarning = warning_on_one_line


def optionParser(test_farm:str, exedir:str):
    """
    Parse command line inputs 
    
    Parse: 
      action to perform
      test/s to run
      executable choice
      number of MPI processes
      number of OMP threads

    :param test_farm: test farm directory
    :param exedir: exciting executable directory 
    :return input_options: Dictionary of parsed command line arguments 
    """
 
    p = ap.ArgumentParser(description=\
        "Usage: python3 runtest.py -a <action> -t <tests> -e <execuatble> -np <NP> -omp <omp> -handle-errors -run-failing-tests")

    help_action = "Defines what action is done. " \
        + "'run' for running tests; " \
        + "'ref' for running references; " \
        + "'clean' for cleaning the test directories; "

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
        + "'exciting_serial' for the serial binary; "\
        + "'exciting_smp' for the shared-memory version; "\
        + "'exciting_mpi' for the binary with MPI parallisation, only; "\
        + "'exciting_mpismp for the binary with MPI amd SMP parallisation;" \
        + "Default is exciting_smp"
    
    p.add_argument('-e',
                   metavar = '--executable',
                   help = help_executable,
                   type = str,
                   default = settings.binary_smp, 
                   choices = settings.binary_names)
    
    p.add_argument('-np',
                   metavar = '--NP',
                   help = "Number of cores for MPI run. Can only be used in " + \
                   "combination with exciting_mpi or exciting_mpismp as executable. " + \
                   "Default is 2 for MPI and MPI+OMP calculations, and 1 for serial or pure OMP", 
                   type = int)
    
    p.add_argument('-omp', 
                   metavar='--ompthreads', 
                   help='Number of OMP threads. ' + \
                   'Default is 2 for OMP and MPI+OMP calculations, and 1 for serial or pure MPI',
                   type = int)

    p.add_argument('-handle-errors',
                    help = 'Allow assertion failures and skips to propagate to the end ' + \
                   'of the test suite. If the option is excluded, the default is to not allow ' + \
                   'error propagation',
                    dest='handle_errors', 
                    default = False, 
                    action='store_true')

    p.add_argument('-run-failing-tests',\
                    help = 'Run tests tagged as failing.' + \
                   'If the option is excluded, the default is not to run failing tests',
                    dest='run_failing_tests', 
                    default = False, 
                    action='store_true')        

    p.add_argument('-make-test', \
                    help='Run tests from Makefile. ' + \
                         'If this option is set, all other options will be ignored ' + \
                         'and the test suite will run all tests with default settings. ' + \
                         'The executable will be chosen from the compiled binaries with the following ' + \
                         'hierarchy: ' + \
                         'exciting_mpismp  >  exciting_smp  >  exciting_mpi  >  exciting_serial .' + \
                         'If excitingtools is not installed, the test suite will provide instructions on how to install the package.',
                    dest='make_test',
                    default=False,
                    action='store_true')


    args = p.parse_args()

    if args.make_test:
        return set_up_make_test(test_farm, exedir)

    input_options = {'action': args.a, 
                     'tests': args.t, 
                     'handle_errors': args.handle_errors, 
                     'run_failing_tests': args.run_failing_tests
                     }

    build_type = args.e if isinstance(args.e, Build_type) else build_type_str_to_enum[args.e]
    input_options['np'] = args.np if args.np is not None else settings.default_np[build_type]
    input_options['omp'] = args.omp if args.omp is not None else settings.default_threads[build_type]

    if args.e == settings.binary_serial:
        if args.np is not None:
            warnings.warn('Using serial exciting, -np will be ignored.')
        if args.omp is not None:
            warnings.warn('Using serial exciting, -omp will be ignored.')
    if args.e == settings.binary_smp:
         if args.np is not None:
            warnings.warn('Using smp exciting, -np will be ignored.')
    if args.e == settings.binary_purempi:
         if args.omp is not None:
            warnings.warn('Using pure mpi exciting, -omp will be ignored.')
            
    if input_options['tests'][0] == 'all':
        input_options['tests'] = next(os.walk(test_farm))[1]
    else:
        input_options['tests'] = input_options['tests']

    for test in input_options['tests']:
        i = input_options['tests'].index(test)
        if test_farm in test:
            input_options['tests'][i] = os.path.basename(os.path.normpath(test))

    test_dirs = next(os.walk(test_farm))[1]
    for test in input_options['tests']:
        if test not in test_dirs:
            raise ValueError('Could not find %s/ in %s/'%(test, test_farm))

    executable_string = os.path.join(exedir, build_type_enum_to_str[build_type])         
    if build_type in [settings.binary_purempi, settings.binary_mpismp]:
        executable_string = \
            'mpirun -np %i %s' % (input_options['np'], executable_string)
    input_options['executable'] = executable_string

    return input_options 

def set_up_make_test(test_farm:str, exedir:str) -> dict:
    """
    Set up test suite for running via make test. The function will check which binaries were compiled
    and set up the executable with the following hierarchy:

        exciting_mpismp > exciting_smp > exciting_purempi > exciting_serial

    If excitingtools is not installed, the function will do that if the user wishes.
    Then the function sets up the default input_options and runs all test cases.
    :param test_farm: test farm directory
    :param exedir: exciting executable directory 
    :return input_options: Dictionary of parsed command line arguments 
    """
    if 'excitingtools' not in sys.modules:
        while(True):
            answer = input('The test suit requires the python package exctingtools. Do you wish to install excitingtools now? (yes/no) \n')
            if answer.lower() in ['yes', 'y']:
                install_excitingtools()
                break
            elif answer.lower() in ['no', 'n']:
                print('You can install excitingtools later via \n\n' + \
                      '    pip install -e %s \n \n'%os.path.join(os.environ('EXCITINGSCRIPTS'), 'exciting_tools') )
                exit()
    
    input_options = {'action': 'run', 
                     'tests': next(os.walk(test_farm))[1], 
                     'handle_errors': False, 
                     'run_failing_tests': False}

    compiled_binaries = []
    for x in next(os.walk('../bin/'))[2]:
        try:
            compiled_binaries.append(build_type_str_to_enum[x])
        except Exception:
            pass
          
    if settings.binary_mpismp in compiled_binaries:
        build_type = settings.binary_mpismp

    elif settings.binary_smp in compiled_binaries:
        build_type = settings.binary_smp

    elif settings.binary_purempi in compiled_binaries:
        build_type = settings.binary_purempi

    elif settings.binary_serial in compiled_binaries:
        build_type = settings.binary_serial
    
    else:
        raise Exception('Could not find any exciting binary. Running the test suite requires compiling ' + \
                        'one of the following exciting binaries:\n\n' + \
                        '    %s  %s  %s  %s. \n\n'%(settings.binary_mpismp, settings.binary_smp,
                                             settings.binary_purempi, settings.binary_serial))
    
    input_options['np'] = settings.default_np[build_type]
    input_options['omp'] = settings.default_threads[build_type]
    
    executable_string = os.path.join(exedir, build_type_enum_to_str[build_type])         
    if build_type in [settings.binary_purempi, settings.binary_mpismp]:
        executable_string = \
            'mpirun -np %i %s' % (input_options['np'], executable_string)
    input_options['executable'] = executable_string
    
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
    
    print(input_options['executable'])
    species_files = next(os.walk(settings.species))[2]
    action = input_options['action']
    executable_command = input_options['executable']
    executable = executable_command.split('/')[-1]
    print(executable)
    skipped_tests = set_skipped_tests(executable, input_options['run_failing_tests']) 

    if action == "run":
        run_tests(settings.test_farm, 
                 settings.main_output,
                 input_options['tests'],
                 settings.run_dir, 
                 settings.ref_dir,
                 settings.input_file,
                 species_files,
                 settings.init_default,
                 executable_command,
                 input_options['np'],
                 input_options['omp'],
                 settings.max_time,
                 skipped_tests,
                 input_options['handle_errors'])
     
    elif action == "ref":    
        while(True):
            answer = input("Are you sure that you want to rerun reference(s)? " + \
                           "This will replace the old reference(s) with reference(s) " + \
                           "created by the current state of the code and should only " + \
                           "be done if all tests succeed with the current version of the code.(yes/nn)\n")
            if answer.lower() in ['yes', 'y']:
                break
            elif answer.lower() in ['no', 'n']:
                exit()
        warnings.warn("Reference will always be run with %s!"%settings.exe_ref)
        run_references(settings.test_farm,
                      settings.main_output,
                      input_options['tests'],
                      settings.ref_dir,
                      os.path.join(settings.exe_dir, build_type_enum_to_str[settings.exe_ref]),
                      settings.ignored_output,
                      settings.max_time)
    

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
        
        clean_tests(settings.test_farm, 
                   input_options['tests'],
                   settings.run_dir,
                   settings.ref_dir,
                   settings.ignored_output)
        exit()
        
    self_test_results = run_self_tests()

    main(settings, input_options)
