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
from typing import List
import re

from tools.termcolor_wrapper import print_color
from tools.runner.test import run_tests
from tools.runner.reference import run_single_reference
from tools.constants import settings
from tools.selftests import runselftests, test_compare
from tools.parsers import install_excitingtools
from tools.utils import Build_type, build_type_str_to_enum, build_type_enum_to_str
from failing_tests import set_skipped_tests


def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return '%s: %s \n' % (category.__name__, message)


warnings.formatwarning = warning_on_one_line


def option_parser(test_farm: str, exe_dir: str):
    """
    Parse command line inputs 
    
    Parse: 
      action to perform
      test/s to run
      executable choice
      number of MPI processes
      number of OMP threads

    :param test_farm: test farm directory
    :param exe_dir: exciting executable directory 
    :return input_options: Dictionary of parsed command line arguments 
    """

    p = ap.ArgumentParser(
        description="Usage: python3 runtest.py -a <action> -t <tests> -e <executable> -np <NP> -omp <omp> -handle-errors -run-failing-tests")

    # TODO(Alex) ref functionality should be moved
    help_action = "Defines what action is done. " \
                  + "'run' for running tests; " \
                  + "'ref' for running references; "

    p.add_argument('-a',
                   metavar='--action',
                   help=help_action,
                   type=str,
                   choices=settings.action_choices,
                   default='run')

    help_test = "Tests to run. This can be some partial string match, such as `PBE`, or a full path " \
                "`groundstate/LDA_PW-PbTiO3`. The test suite will run any test which contain this " \
                "string in its name. Typing 'all' will run all test cases, and is the default behaviour."

    p.add_argument('-t',
                   metavar='--tests',
                   help=help_test,
                   type=str,
                   default=['all'],
                   nargs='+')

    help_executable = "exciting executables. " \
                      + "'exciting_serial' for the serial binary; " \
                      + "'exciting_smp' for the shared-memory version; " \
                      + "'exciting_mpi' for the binary with MPI parallisation, only; " \
                      + "'exciting_mpismp for the binary with MPI amd SMP parallisation;" \
                      + "Default is exciting_smp"

    p.add_argument('-e',
                   metavar='--executable',
                   help=help_executable,
                   type=str,
                   default=settings.binary_smp,
                   choices=settings.binary_names)

    p.add_argument('-np',
                   metavar='--NP',
                   help="Number of cores for MPI run. Can only be used in " +
                        "combination with exciting_mpi or exciting_mpismp as executable. " +
                        "Default is 2 for MPI and MPI+OMP calculations, and 1 for serial or pure OMP",
                   type=int)

    p.add_argument('-omp',
                   metavar='--ompthreads',
                   help='Number of OMP threads. ' +
                        'Default is 2 for OMP and MPI+OMP calculations, and 1 for serial or pure MPI',
                   type=int)

    p.add_argument('-handle-errors',
                   help='Allow assertion failures and skips to propagate to the end ' +
                        'of the test suite. If the option is excluded, the default is to not allow ' +
                        'error propagation',
                   dest='handle_errors',
                   default=False,
                   action='store_true')

    p.add_argument('-run-failing-tests',
                   help='Run tests tagged as failing.' +
                        'If the option is excluded, the default is not to run failing tests',
                   dest='run_failing_tests',
                   default=False,
                   action='store_true')

    p.add_argument('-make-test',
                   help='Run tests from Makefile. ' +
                        'If this option is set, all other options will be ignored ' +
                        'and the test suite will run all tests with default settings. ' +
                        'The executable will be chosen from the compiled binaries with the following ' +
                        'hierarchy: ' +
                        'exciting_mpismp  >  exciting_smp  >  exciting_mpi  >  exciting_serial .' +
                        'If excitingtools is not installed, the test suite will provide instructions on how to install the package.',
                   dest='make_test',
                   default=False,
                   action='store_true')

    args = p.parse_args()

    if args.make_test:
        return set_up_make_test(test_farm, exe_dir)

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

    all_test_dirs = get_test_directories(test_farm)

    if input_options['tests'][0] == 'all':
        input_options['tests'] = all_test_dirs
    else:
        tests_to_run = partial_test_name_matches(all_test_dirs, input_options['tests'])
        if not tests_to_run:
            raise ValueError('Could not find any full or partial string matches to %s for test names in '
                             '%s subdirectories' % (input_options['tests'], test_farm))
        input_options['tests'] = tests_to_run

    executable_string = os.path.join(exe_dir, build_type_enum_to_str[build_type])

    if not os.path.isfile(executable_string):
        quit('exciting binary not found in' + exe_dir)

    if build_type in [settings.binary_purempi, settings.binary_mpismp]:
        executable_string = 'mpirun -np %i %s' % (input_options['np'], executable_string)
    input_options['executable'] = executable_string

    return input_options


# TODO(Alex) Move in a new MR
def get_test_directories(test_farm_root: str, basename=False) -> List[str]:
    """
    Get test directories from the test farm root.

    Test farm has the directory structure:
     test_farm/method/test_directory

    :param str test_farm_root: Test farm root directory name
    :param bool basename: Only return test directory base names
    :return List[str] test_dirs: List of test directories, given relative to test_farm_root
    if basename is false.
    """

    method_dirs = next(os.walk(test_farm_root))[1]

    if basename:
        test_dirs = []
        for method_dir in method_dirs:
            method_dir = os.path.join(test_farm_root, method_dir)
            test_dirs += next(os.walk(method_dir))[1]
        return test_dirs

    else:
        full_dirs = []
        for method_dir in method_dirs:
            method_dir = os.path.join(test_farm_root, method_dir)
            test_dirs = next(os.walk(method_dir))[1]
            full_dirs += [os.path.join(method_dir, t) for t in test_dirs]
        return full_dirs


def partial_test_name_matches(all_tests: List[str], input_tests: List[str]) -> List[str]:
    """
    Given a list of strings, return any full or partial matches in the test_farm
    subdirectories.

    :param List[str] all_tests: All tests in the test_farm
    :param List[str] input_tests: List of partial test names, provided from input
    :return List[str] tests_to_run: List of tests to run
    """

    all_tests_str = "\n".join(test for test in all_tests)

    tests_to_run = []
    for test in input_tests:
        tests_to_run += re.findall("^.*" + test + ".*$", all_tests_str, re.MULTILINE)

    return tests_to_run


def set_up_make_test(test_farm: str, exe_dir: str) -> dict:
    """
    Set up test suite for running via make test. The function will check which binaries were compiled
    and set up the executable with the following hierarchy:

        exciting_mpismp > exciting_smp > exciting_purempi > exciting_serial

    If excitingtools is not installed, the function will do that if the user wishes.
    Then the function sets up the default input_options and runs all test cases.

    :param str test_farm: test farm directory
    :param str exe_dir: exciting executable directory
    :return dict input_options: Dictionary of parsed command line arguments
    """
    if 'excitingtools' not in sys.modules:
        while True:
            answer = input(
                'The test suit requires the python package exctingtools. Do you wish to install excitingtools now? (yes/no) \n')
            if answer.lower() in ['yes', 'y']:
                install_excitingtools()
                break
            elif answer.lower() in ['no', 'n']:
                print('You can install excitingtools later via \n\n' + \
                      '    pip install -e %s \n \n' % os.path.join(os.environ('EXCITINGSCRIPTS'), 'exciting_tools'))
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
        raise Exception('Could not find any exciting binary. Running the test suite requires compiling ' 
                        'one of the following exciting binaries:\n\n'
                        '    %s  %s  %s  %s. \n\n' % (settings.binary_mpismp, settings.binary_smp,
                                                      settings.binary_purempi, settings.binary_serial))

    input_options['np'] = settings.default_np[build_type]
    input_options['omp'] = settings.default_threads[build_type]

    executable_string = os.path.join(exe_dir, build_type_enum_to_str[build_type])
    if build_type in [settings.binary_purempi, settings.binary_mpismp]:
        executable_string = 'mpirun -np %i %s' % (input_options['np'], executable_string)
    input_options['executable'] = executable_string

    return input_options


def main(settings: namedtuple, input_options: dict):
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

    skipped_tests = set_skipped_tests(executable, input_options['run_failing_tests'])

    if action == "run":
        run_tests(settings.main_output,
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

    # TODO(Alex) This functionality should be moved 
    elif action == "ref":
        while True:
            answer = input("Are you sure that you want to rerun reference(s)? " +
                           "This will replace the old reference(s) with reference(s) " +
                           "created by the current state of the code and should only " +
                           "be done if all tests succeed with the current version of the code.(yes/nn)\n")
            if answer.lower() in ['yes', 'y']:
                break
            elif answer.lower() in ['no', 'n']:
                exit()

        warnings.warn("Reference will always be run with %s!" % settings.exe_ref)

        executable = os.path.join(settings.exe_dir, build_type_enum_to_str[settings.exe_ref])

        for test in input_options['tests']:
            reference_dir = os.path.join(test, settings.ref_dir)
            run_single_reference(reference_dir, executable, settings)


def run_self_tests(verbosity=3):
    """
    Test the application test suite functions 
    Execute self tests before the test suite
   
    :return: result  
    """
    print('Run self tests:')
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    suite.addTests(loader.loadTestsFromModule(runselftests))
    suite.addTests(loader.loadTestsFromModule(test_compare))

    runner = unittest.TextTestRunner(verbosity=verbosity)
    result = runner.run(suite)

    if len(result.errors) == 0:
        print_color('Self tests SUCCESS', 'green')
    else:
        print_color('Self tests FAIL', 'red')
        quit('Test suite has quit')

    return result


if __name__ == "__main__":
    input_options = option_parser(settings.test_farm, settings.exe_dir)
    self_test_results = run_self_tests()
    main(settings, input_options)
