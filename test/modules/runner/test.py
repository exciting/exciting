import os
from typing import List, Union, Tuple
import sys
import json
from copy import deepcopy

from excitingtools.dict_utils import delete_nested_key

from ..constants import files_under_test, settings, methods_moved_to_json, keys_to_remove
from ..parsers import ErroneousFileError, parser_chooser
from ..infrastructure import create_run_dir, copy_exciting_input, get_test_from_init, flatten_directory

from ..termcolor_wrapper import print_color
from ..tester.compare import ErrorFinder
from ..tester.test import from_init
from ..tester.report import Report, test_suite_summary, skipped_test_summary, timing_summary, TestResults, \
    SummariseTests
from ..tester.failure import Failure, Failure_code

from .execute import execute


def read_output_file(file_name: str) -> Union[dict, Failure]:
    """
    Read an exciting output file

    If the subdirectory in file_name is not ref

    :param str file_name: File name
    :return dict ref_data: Reference data
    """
    sub_dir = file_name.split('/')[-2]
    assert sub_dir in ['ref', 'run'], "read_output_file: Subdirectory in which file exists is not 'ref' or 'run' "

    try:
        data = parser_chooser(file_name)
        remove_untested_keys(data, file_name, keys_to_remove)
        return data
    except OSError:
        failure_code = {'ref': Failure_code.REFERENCE, 'run': Failure_code.RUN}
        return Failure(test_name=file_name, failure_code=failure_code[sub_dir])
    except ErroneousFileError:
        return Failure(test_name=file_name, failure_code=Failure_code.ERRORFILE)


def remove_untested_keys(data: dict, file_name: str, keys_to_remove: dict) -> dict:
    """
    Remove keys not to be tested from data dictionary

    :param dict data: Test or reference data
    :param str file_name: File name
    :param dict keys_to_remove: dictionary of keys not to be tested

    """
    
    # Remove all scf loops except for the last one
    if 'INFO.OUT' in file_name:
        scl_indices = [int(item) for item in list(data['scl'].keys())]
        last_scl = max(scl_indices)
        data['scl'] = data['scl'][str(last_scl)]

    #remove untested keys
    for key_chain in keys_to_remove[file_name.split('/')[-1]]:
        delete_nested_key(data, key_chain)

        
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

        ref_data = read_output_file(reference_file)
        test_data = read_output_file(test_file)

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


# TODO(Alex) Issue 99.  Function to remove once test suite migrates to JSON
def compare_outputs_xml(output_files: dict, test_dir: str, ref_dir: str, run_dir: str, report: Report):
    """
    Compare calculation outputs (output_files) to reference data, with tolerances defined in init.xml

    :param dict output_files: File names and associated tolerances to validate for a given method, specified in init.xml
    :param str test_dir: Test directory
    :param str ref_dir: Reference directory (relative to test_dir)
    :param str run_dir: Run directory (relative to test_dir)
    :param Report report: Report object, summarising test failures
    """
    # Loop over all files to test and compare them to their references
    for test in output_files:
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
    return report


def run_single_test_json(main_out: str, test_dir: str, run_dir: str, ref_dir: str, input_file: str,
                         species_files: List[str], executable: str, max_time: int, handle_errors: bool) -> TestResults:
    """
    Runs a test case and compares output files under test against references.

    :param main_out:        main output file of the exciting calculation
    :param test_dir:        test case that will be run
    :param run_dir:         name of the run directory of the test case
    :param ref_dir:         name of the ref directory of the test case
    :param input_file:      exciting input file
    :param species_files:   list of species files
    :param executable:      executable command
    :param max_time:        max time a job is allowed to run for before being killed
    :param handle_errors:   Whether or not failures and skips are allowed to propagate

    :return TestResults test_results: Test case results object.
    """
    head_tail = os.path.split(test_dir)
    test_name = head_tail[1]
    method = head_tail[0].split('/')[-1]
    full_ref_dir = os.path.join(test_dir, ref_dir)
    full_run_dir = os.path.join(test_dir, run_dir)

    print('Run test %s:' % test_name)

    json_tolerances = load_tolerances(full_ref_dir)
    # TODO(A/B/H) Issue 100. Update ErrorFinder class to use tol AND units in data comparison
    #  As an alternative to stripping the units (could pass the comparison function as an argument)
    # See function documentation for a description of the issue
    just_tolerances = strip_tolerance_units(json_tolerances)

    output_files = files_under_test[method]
    create_run_dir(test_dir, run_dir)
    copy_exciting_input(full_ref_dir, full_run_dir, species_files, input_file)

    run_success, err_mess, timing = execute(full_run_dir, executable, main_out, max_time)
    test_results = TestResults(test_dir, run_success, timing)

    flatten_directory(full_run_dir)

    if run_success:
        test_results_dict = compare_outputs_json(full_run_dir, full_ref_dir, output_files, just_tolerances)
        test_results.set_results(test_results_dict)

    test_results.print_results()
    test_results.assert_errors(handle_errors)

    return test_results


# TODO(Alex) Issue 99.  Function to remove once test suite migrates to JSON
def run_single_test_xml(main_out: str, test_dir: str, run_dir: str, ref_dir: str, input_file: str,
                        species_files: List[str], init_default: str, executable: str, max_time: str, timing: dict,
                        handle_errors: bool) -> Report:
    """
    Runs a single test.

    :param main_out:        main output file of the exciting calculation
    :param test_dir:        test case that will be run
    :param run_dir:         name of the run directory of the test case
    :param ref_dir:         name of the ref directory of the test case
    :param input_file:      exciting input file
    :param species_files:   list of species files
    :param init_default:    location of the default init.xml
    :param executable:      executable command
    :param timing:          test run times in seconds
    :param handle_errors:   Whether or not failures and skips are allowed to propagate

    :return Report report:  Report instance of the test
    """
    test_name = os.path.basename(test_dir)
    print('Run test %s:' % test_name)

    full_ref_dir = os.path.join(test_dir, ref_dir)
    full_run_dir = os.path.join(test_dir, run_dir)

    # Parse tolerances and output_files
    init = get_test_from_init(test_dir, init_default)
    name = init['name']
    description = init['description']
    output_files = init['tests']

    report = Report(name, description)

    create_run_dir(test_dir, run_dir)

    copy_exciting_input(full_ref_dir, full_run_dir, species_files, input_file)

    # Run exciting.
    # If the run fails, then an error message is added and all other init['tests'] are skipped.
    run_success, err_mess, timing[init['name']] = execute(full_run_dir, executable, main_out, max_time)

    if not run_success:
        report.runFailed(test_name, err_mess)
        print('Time (s): %.1f' % timing[name])
        report.writeToTerminal()
        report.assert_errors(handle_errors)
        return report
    print('Run succeeded')
    print(full_run_dir)

    flatten_directory(full_run_dir)

    report = compare_outputs_xml(output_files, test_dir, ref_dir, run_dir, report)
    print('Time (s): %.1f' % timing[name])
    report.writeToTerminal()
    report.assert_errors(handle_errors)

    return report


def split_test_list_according_to_tol_format(test_list: List[str]) -> tuple:
    """
    Split the test list according to which methods use XML tolerances and which
    use JSON tolerances.

    Expect test_list to contain elements of the form: test_farm/groundstate/LDA_PW-PbTiO3
    however, it does not affect the routine. i.e.

    `if 'groundstate' in "test_farm/groundstate"`

    will give the same result as

    `if 'groundstate' in "groundstate"`

    TODO(Alex) Issue 99.  Function to remove once test suite migrates to JSON
         Gradually extend HARDCODED definition of which test methods are run with JSON
         Once all methods are performed with JSON, remove the XML calls and delete this routine.

    :param List[str] test_list: List test cases that will be run.
    :return tuple test_list_json, test_list_xml: List of test cases for JSON and XML, respectively.
    """
    test_list_xml = []
    test_list_json = []

    def test_matches_method(test_method: str, methods: List[str]) -> bool:
        for m in methods:
            if test_method == m:
                return True
        return False

    for full_test_name in test_list:
        head_tail = os.path.split(full_test_name)
        test_name = head_tail[1]
        test_method = head_tail[0].split('/')[-1]

        if test_matches_method(test_method, methods_moved_to_json):
            test_list_json.append(full_test_name)
        else:
            test_list_xml.append(full_test_name)

    return test_list_json, test_list_xml


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

    test_list, removed_tests = remove_tests_to_skip(test_list, skipped_tests)
    # TODO(Alex) Issue 99. Remove once all method tolerances are moved to JSON
    test_list_json, test_list_xml = split_test_list_according_to_tol_format(test_list)

    # JSON
    print_color('\n\nTests using JSON-based tolerances', 'blue')
    report = SummariseTests()
    for test_name in test_list_json:
        test_results = run_single_test_json(main_out, test_name, run_dir, ref_dir, input_file, species_files,
                                            executable, max_time, handle_errors)
        report.add(test_results)

    report.print()
    skipped_test_summary(removed_tests)
    report.print_timing()
    assert report.n_failed_test_cases == 0, "Some test suite cases failed."

    # XML-based test suite
    print_color('\n\nTests using XML-based tolerances', 'blue')
    timing = {}
    test_suite_report = []
    for test_name in test_list_xml:
        test_suite_report.append(
            run_single_test_xml(main_out,
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
    # TODO(Alex). Issue 98. Move Printing of Skipped Test Cases to SummariseTests class
    skipped_test_summary(removed_tests)
    timing_summary(timing, verbose=True)
    assert all_asserts_succeeded, "Some test suite assertions failed or were skipped over"

    return


def remove_tests_to_skip(test_list: List[str], skipped_tests: List[dict]) -> Tuple[List[str], List[dict]]:
    """
    Remove tests given in 'skipped_tests' from the test suite, for a specific executable choice.

    This is useful if a particular test crashes or hangs, and needs to be
    debugged BUT shouldn't cause the test suite to report a failure.

    Notes
    ------------
    Using sets is probably faster but they won't preserve ordering. Could use ordered sets
    but would need to support that package

    Arguments
    -------------
    :param List[str] test_list: List of test names to run.
    :param List[dict] skipped_tests: List of tests from the test suite that are marked as "to skip",
    in failing_tests.py

    :return List[str] tests_to_run: Some list of tests to run
    :return List[dict] removed_tests: List of tests removed from test_list, where each entry is a dict
    (like with skipped_tests).
    """
    test_farm_root = settings.test_farm

    tests_to_skip = [os.path.join(test_farm_root, test['name']) for test in skipped_tests]
    comment = {os.path.join(test_farm_root, test['name']): test['comment'] for test in skipped_tests}

    tests_to_run = []
    removed_tests = []

    for test in test_list:
        if test not in tests_to_skip:
            tests_to_run.append(test)
        else:
            removed_tests.append({'name': test, 'comment': comment[test]})

    return tests_to_run, removed_tests
