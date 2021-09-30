import os
from typing import List, Union
import sys
import json
from copy import deepcopy

from ..parsers import ErroneousFileError, parser_chooser
from ..infrastructure import create_run_dir, copy_exciting_input, get_test_from_init, flatten_directory

from ..tester.compare import ErrorFinder
from ..tester.test import from_init
from ..tester.report import Report, test_suite_summary, skipped_test_summary, timing_summary
from ..tester.failure import Failure, Failure_code

from .execute import execute


def read_output_file(file_name: str, calc_type: str) -> Union[dict, Failure]:
    """
    Read an exciting output file

    :param str file_name: File name
    :param str calc_type: Calculation type (such that error handling is correct)
     Only valid strings are 'ref' and 'run' (this could be implemented better)
    :return dict ref_data: Reference data
    """
    assert calc_type in ['ref', 'run'], 'Specify "ref" or "run" when calling "read_output_file"'

    try:
        ref_data = parser_chooser(file_name)
        return ref_data
    except OSError:
        failure_code = {'ref': Failure_code.REFERENCE, 'run': Failure_code.RUN}
        return Failure(test_name=file_name, failure_code=failure_code[calc_type])
    except ErroneousFileError:
        return Failure(test_name=file_name, failure_code=Failure_code.ERRORFILE)


def load_tolerances(directory: str) -> dict:
    """
    Get the tolerances from any files in 'directory' of the form `*tol*.json`

    :param str directory: Directory containing `*tol*.json` file
    :return dict tolerances: Dictionary of tolerances
    """

    file_names = next(os.walk(directory))[2]
    tolerance_files = [f for f in file_names if ('tol' in f) and ('.json' in f)]

    if not tolerance_files:
        sys.exit('No tolerance files in ' + directory)

    tolerances = {}
    for file in tolerance_files:
        with open(os.path.join(directory, file)) as fid:
            tolerances.update(json.load(fid))

    return tolerances


def strip_tolerance_units(json_tolerance: dict) -> dict:
    """
    Strip the units from the json_tolerances:

    json_tolerance[file_name][key] =  {'tol': 1e-08, 'unit': 'Bohr', 'comment': 'text'}
    to
    just_tolerance[file_name][key] =  1e-08
    """

    just_tolerances = {}
    for file_name, tolerances in json_tolerance.items():
        tmp = {}
        for key, entry in tolerances.items():
            tmp[key] = entry['tol']
        just_tolerances[file_name] = deepcopy(tmp)
        del tmp
    return just_tolerances


def compare_outputs_json(run_dir: str, ref_dir: str, output_files: List[str], tolerance: dict) -> dict:
    """
    For a given test case, compare each file specified `output_files` to reference data.

    TODO(A/B/H) Issue 97. Unevaluated Test Data Keys
     Check for missing keys in tolerances, from what's returned by ErrorFinder

    :param List[str] output_files: Files to validate for a given method, specified in init.xml
    :param str ref_dir: Reference directory
    :param str run_dir: Run directory
    :param List[str] output_files: List of files to compare against for this calculation
    :param dict tolerance: Tolerances, of form {'method': {'key1': value1, ...}}

    :return dict test_results: Dictionary of errors per file, for the test case.
    """
    test_results = {}
    for file in output_files:

        reference_file = os.path.join(ref_dir, file)
        test_file = os.path.join(run_dir, file)

        ref_data = read_output_file(reference_file, 'ref')
        test_data = read_output_file(test_file, 'run')

        # Handle IO errors
        if isinstance(ref_data, Failure):
            test_results[file] = ref_data
            continue
        if isinstance(test_data, Failure):
            test_results[file] = test_data
            continue

        try:
            file_tolerances = tolerance[file]
        except KeyError:
            raise KeyError("File name does not match the file_name_key in tolerance file."
                           "Most likely the key defined in tolerance templates is not equal"
                           "to the output file name")

        file_results = ErrorFinder(test_data, ref_data, file_tolerances, file_name=file)
        test_results[file] = deepcopy(file_results)

    return test_results


############################################

def run_single_test(main_out: str, test_dir: str, run_dir: str, ref_dir: str, input_file: str,
                    species_files: List[str], init_default: str, executable: str, max_time: str, timing: dict,
                    handle_errors: bool):
    """
    Runs a single test.
    :param main_out:        main output file of the exciting calculation
    :param test_dir:       test case that will be run
    :param run_dir:         name of the run directory of the test case
    :param ref_dir:         name of the ref directory of the test case
    :param input_file:      exciting input file
    :param species_files:   list of species files
    :param init_default:    location of the default init.xml
    :param executable:     executable command
    :param timing:          test run times in seconds
    :param handle_errors:   Whether or not failures and skips are allowed to propagate

    Output:
        report          object      report instance of the test
    """
    test_name = os.path.basename(test_dir)
    print('Run test %s:' % test_name)

    # parse init file --> parse tolerances and tests
    init = get_test_from_init(test_dir, init_default)
    name = init['name']
    description = init['description']
    tests = init['tests']

    report = Report(name, description)

    create_run_dir(test_dir, run_dir)

    copy_exciting_input(os.path.join(test_dir, ref_dir),
                        os.path.join(test_dir, run_dir),
                        species_files,
                        input_file)

    # Run exciting.
    # If the run fails, then an error message is added
    # and all other init['tests'] are skipped.
    run_success, err_mess, timing[init['name']] = execute(os.path.join(test_dir, run_dir),
                                                          executable,
                                                          main_out,
                                                          max_time)
    if not run_success:
        report.runFailed(test_name, err_mess)
        print('Time (s): %.1f' % timing[name])
        report.writeToTerminal()
        report.assert_errors(handle_errors)
        return report
    print('Run succeeded')
    print(os.path.join(test_dir, run_dir))

    # Flatten run directory
    flatten_directory(os.path.join(test_dir, run_dir))

    # Loop over all files to test and compare them to their reference
    for test in tests:
        file_name = test['file']
        test_results = from_init(test)

        # read reference data
        ref_path = os.path.join(test_dir, ref_dir, file_name)
        try:
            ref_data = parser_chooser(ref_path)
        except OSError:
            test_results.append(Failure(Failure_code.REFERENCE, err_msg=file_name))
            report.collectTest(test_results)
            continue

        # read run data
        run_path = os.path.join(test_dir, run_dir, file_name)
        try:
            run_data = parser_chooser(run_path)
        except OSError:
            test_results.append(Failure(Failure_code.FILENOTEXIST, err_msg=file_name))
            report.collectTest(test_results)
            continue
        except ErroneousFileError:
            test_results.append(Failure(Failure_code.ERRORFILE, err_msg=file_name))
            report.collectTest(test_results)
            continue

        # compare reference data to run data
        test_results.evaluate(ref_data, run_data, test_dir)

        report.collectTest(test_results)

    print('Time (s): %.1f' % timing[name])
    report.writeToTerminal()
    report.assert_errors(handle_errors)

    return report


def run_tests(main_out: str, test_list: List[str], run_dir: str, ref_dir: str,
              input_file: str, species_files: List[str], init_default: str,
              executable: str, np: int, omp: int, max_time: int,
              skipped_tests: List[str], handle_errors: bool):
    """
    Runs tests in test_list (see run_single_test).
    :param main_out:          main output file of the exciting calculation
    :param test_list:         list of string test cases that will be run
    :param run_dir:           name of the run directory of the test case
    :param ref_dir:           name of the ref directory of the test case
    :param input_file:        exciting input file
    :param species_files:     list of species files
    :param init_default:      location of the default init.xml
    :param executable:        executable command for the exciting run
    :param np:                number of MPI processes
    :param omp:               number of OMP threads
    :param max_time:          max time before a job is killed
    :param skipped_tests:     list of tests to skip
    :param handle_errors:     Whether or not failures and skips are allowed to propagate
    """
    if 'exciting_serial' in executable:
        print('Run tests with exciting_serial.')
    elif 'exciting_smp' in executable:
        print('Run tests with exciting_smp with %i open MP threads.' % (omp))
    elif 'exciting_mpismp' in executable:
        print('Run tests with exciting_mpismp with %i open MP threads and %i MPI processes.' % (omp, np))
    elif 'exciting_purempi' in executable:
        print('Run tests with exciting_purempi %i MPI processes.' % np)

    test_list = remove_tests_to_skip(test_list, skipped_tests)
    timing = {}
    test_suite_report = []

    for test_name in test_list:
        test_suite_report.append(
            run_single_test(main_out,
                            test_name,
                            run_dir,
                            ref_dir,
                            input_file,
                            species_files,
                            init_default,
                            executable,
                            max_time,
                            timing,
                            handle_errors)
        )

    all_asserts_succeeded = test_suite_summary(test_suite_report)
    skipped_test_summary(skipped_tests)
    timing_summary(timing, verbose=True)
    assert all_asserts_succeeded, "Some test suite assertions failed or were skipped over"

    return


def remove_tests_to_skip(all_tests: List[str], skipped_tests: List[str]) -> List[str]:
    """
    Remove tests given in 'skipped_tests' from the test suite,
    for a specific executable choice.

    This is useful if a particular test crashes or hangs, and needs to be
    debugged BUT shouldn't cause the test suite to report a failure.

    :param all_tests:     list of all test names
    :param skipped_tests: list of tests to skip. Each entry is a dict

    :return tests_to_run: list of tests to run (with skipped_tests removed)
    """

    tests_to_skip = [test['name'] for test in skipped_tests]
    # Using sets is probably faster but they won't preserve ordering. Could use ordered sets
    tests_to_run = []
    for test in all_tests:
        if os.path.basename(test) not in tests_to_skip:
            tests_to_run.append(test)

    return tests_to_run
