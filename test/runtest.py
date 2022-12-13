"""
Parser command-line arguments and run the test suite.

To run all tests in smp, type `python3 runtest.py`
To run all tests in mpi and smp, type `python3 runtest.py -e exciting_mpismp`
For more details, type `python3 runtest.py --help`
"""
import argparse as ap
import warnings
import os
from typing import List, Union

from src.exciting_settings.constants import settings as exciting_settings, Defaults
from src.runner.execute import set_job_environment
from src.runner.profile import BuildType, build_type_str_to_enum, build_type_enum_to_str
from src.runner.set_tests import get_test_directories, partial_test_name_matches
from src.runner.test import run_tests
from src.tester.report import skipped_test_summary
from src.runner.set_tests import get_test_cases_from_config, set_tests_to_run


def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return '%s: %s \n' % (category.__name__, message)


warnings.formatwarning = warning_on_one_line


def option_parser(settings: Defaults):
    """
    Parse command line inputs

    :param Defaults settings: exciting defaults.
    :return input_options: Dictionary of parsed command line arguments.
    """

    p = ap.ArgumentParser(
        description="Usage: python3 runtest.py -t <tests> -e <executable> -np <NP> -omp <omp> "
                    "-handle-errors -run-failing-tests -repeat-tests <N>")

    help_action = "Defines what action is done. " \
                  + "'run' for running tests (default); " \
                  + "'ref' for (re)running all references; "

    help_test = "Tests to run. This can be some partial string match, such as `PBE`, or a full path " \
                "`groundstate/LDA_PW-PbTiO3`. The test suite will run any test which contain this " \
                "string in its name. Typing 'all' will run all test cases, and is the default behaviour."

    p.add_argument('-t',
                   metavar='--tests',
                   help=help_test,
                   type=str,
                   default=[],
                   nargs='+')

    help_executable = "exciting executables. " \
                      + "'exciting_serial' for the serial binary; " \
                      + "'exciting_smp' for the shared-memory version; " \
                      + "'exciting_purempi' for the binary with MPI parallelisation, only; " \
                      + "'exciting_mpismp for the binary with MPI amd SMP parallelisation;" \
                      + "Default is exciting_smp"

    p.add_argument('-e',
                   metavar='--executable',
                   help=help_executable,
                   type=str,
                   default=settings.binary_smp,
                   choices=settings.binary_names)

    p.add_argument('-np',
                   metavar='--NP',
                   help="Number of cores for MPI run. Can only be used in combination with exciting_mpi or "
                        "exciting_mpismp as executable. Default is 2 for MPI and MPI+OMP calculations, "
                        "and 1 for serial or pure OMP",
                   type=int)

    p.add_argument('-omp',
                   metavar='--ompthreads',
                   help='Number of OMP threads. ' +
                        'Default is 2 for OMP and MPI+OMP calculations, and 1 for serial or pure MPI',
                   type=int)

    p.add_argument('-handle-errors',
                   help='Allow assertion failures to propagate to the end of the test suite. '
                        'If the option is excluded, the default is to not allow error propagation',
                   dest='handle_errors',
                   default=False,
                   action='store_true')

    p.add_argument('-run-failing-tests',
                   help='Run tests tagged as failing.' +
                        'If the option is excluded, the default is not to run failing tests',
                   dest='run_failing_tests',
                   default=False,
                   action='store_true')

    p.add_argument('-repeat-tests',
                   dest='repeat_tests',
                   help='Number of times to repeat any test specified in failing_tests.py/repeat_tests list.'
                        'Primarily intended for flakey tests running in the CI',
                   type=int,
                   default=0)

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

    input_options = {'tests': set_test_names_from_cmd_line(settings.test_farm, args.t),
                     'handle_errors': args.handle_errors,
                     'run_failing_tests': args.run_failing_tests,
                     'repeat_tests': args.repeat_tests
                     }

    if args.make_test:
        return set_up_make_test(settings, input_options)

    build_type = set_build_type(args, settings)
    input_options['np'] = args.np if args.np is not None else settings.default_np[build_type]
    input_options['omp'] = str(args.omp) if args.omp is not None else str(settings.default_threads[build_type])
    input_options['executable'] = set_execution_str(build_type, input_options['np'], settings)
    input_options['mkl_threads'] = set_mkl_threads_from_env()

    return input_options


def set_build_type(cmd_line_args, settings: Defaults) -> BuildType:
    """
    Determine the build type, given a command line argument string.

    Also checks MPI and OMP processes requested, given the BuildType.

    :param cmd_line_args: Command line arguments.
    :param Defaults settings: Default exciting settings. Used to obtain valid build types
    :return BuildType build_type: Build type enum.
    """
    build_type = cmd_line_args.e if isinstance(cmd_line_args.e, BuildType) else build_type_str_to_enum[cmd_line_args.e]

    if cmd_line_args.e == settings.binary_serial:
        if cmd_line_args.np is not None:
            warnings.warn('Using serial exciting, -np will be ignored.')
        if cmd_line_args.omp is not None:
            warnings.warn('Using serial exciting, -omp will be ignored.')

    if cmd_line_args.e == settings.binary_smp:
        if cmd_line_args.np is not None:
            warnings.warn('Using smp exciting, -np will be ignored.')

    if cmd_line_args.e == settings.binary_purempi:
        if cmd_line_args.omp is not None:
            warnings.warn('Using pure mpi exciting, -omp will be ignored.')

    return build_type


def set_test_names_from_cmd_line(test_farm: str, input_tests: List[str]) -> List[str]:
    """
    Set test names to run from specified command line arguments.

    If nothing is specified on the command line, an empty list is returned.
    If partial test names or directories are specified, any string matches are returned.

    :param str test_farm: Relative path to test farm directory, from the test/ root.
    :param List[str] input_tests: Test names or partial string matches to test names.
    :return List[str] tests_to_run: Test names, prepended by path w.r.t. test/ directory.
    """
    if not input_tests:
        return []

    all_test_dirs = get_test_directories(test_farm)
    tests_to_run = partial_test_name_matches(all_test_dirs, input_tests)

    if not tests_to_run:
        raise ValueError(f'Could not find any full or partial string matches to {input_tests} for test names in '
                         f'{test_farm} subdirectories')
    return tests_to_run


def set_mpi_command() -> str:
    """ Set the MPI command.

    mpirun (at least for openMPI) cannot automatically run as root.
    Docker executes commands as root, so one needs to append a flag to mpirun.

    :return mpi_cmd: mpi run command.
    """
    mpi_cmd = "mpirun"
    if os.getenv('OPENMPI_IN_DOCKER') is not None:
        mpi_cmd += " --allow-run-as-root"
    return mpi_cmd


def set_execution_str(build_type: BuildType, np: int, settings: Defaults) -> str:
    """
    Set the execution string.

    :param BuildType build_type: Build type of program binary
    :param int np: Number of MPI processes.
    :param Defaults settings: Default exciting settings. Used to obtain valid build types and executable location.
    :return str executable_string: Execution string.
    """
    executable_string = os.path.join(settings.exe_dir, build_type_enum_to_str[build_type])

    if not os.path.isfile(executable_string):
        raise FileNotFoundError(f'Could not find an exciting binary in {settings.exe_dir}')

    if build_type in [settings.binary_purempi, settings.binary_mpismp]:
        mpi_cmd = set_mpi_command()
        executable_string = f'{mpi_cmd} -np {np} {executable_string}'

    return executable_string


def set_mkl_threads_from_env() -> Union[str, None]:
    """ Set number of threads to be used by MKL from the environment.

    :return mkl_threads: Number of threads to be used by MKL
    """
    mkl_threads = os.getenv("MKL_NUM_THREADS")
    if mkl_threads is not None:
        return str(mkl_threads)
    return mkl_threads


def set_up_make_test(settings: Defaults, input_options: dict) -> dict:
    """
    Set the test suite options for running via `make test`.

    The function will check which binaries are compiled and select the executable according to the following hierarchy:

        exciting_mpismp > exciting_smp > exciting_purempi > exciting_serial

    then overwrite the appropriate execution defaults.

    :param Defaults settings: Default exciting settings.
    :param dict input_options: Dictionary of parsed command line arguments.
    :return dict input_options: Inputs, with appropriate defaults overwritten according to the binary selected.
    """
    compiled_binaries = []
    for x in next(os.walk('../bin/'))[2]:
        try:
            compiled_binaries.append(build_type_str_to_enum[x])
        except KeyError:
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
        raise FileNotFoundError('Could not find any exciting binary. Running the test suite requires compiling '
                                'one of the following exciting binaries:\n\n'
                                '    %s  %s  %s  %s. \n\n' % (settings.binary_mpismp, settings.binary_smp,
                                                              settings.binary_purempi, settings.binary_serial))

    input_options['np'] = settings.default_np[build_type]
    input_options['omp'] = settings.default_threads[build_type]
    input_options['executable'] = set_execution_str(build_type, input_options['np'], settings)
    input_options['mkl_threads'] = set_mkl_threads_from_env()

    return input_options


def main(settings: Defaults, input_options: dict):
    """
    Run test suite

    :param settings: Default settings for environment variables, paths
                     and run-time parameters. 
    :param input_options: parsed command-line arguments
    """
    print('Run test suite:')
    print('Binary: ', input_options['executable'])

    my_env = set_job_environment({key: input_options[key] for key in ['omp', 'mkl_threads']})
    tests_in_config = get_test_cases_from_config(settings)
    tests = set_tests_to_run(settings, input_options, tests_in_config)
    report = run_tests(settings, input_options, tests, my_env)

    report.print()
    skipped_test_summary(tests.failing)
    print("Note, test names in groups that are not run are not printed.")
    report.print_timing()
    assert report.n_failed_test_cases == 0, "Some test suite cases failed."


if __name__ == "__main__":
    input_options = option_parser(exciting_settings)
    main(exciting_settings, input_options)
