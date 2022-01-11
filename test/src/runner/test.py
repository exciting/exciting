"""
Module to execute test cases
"""
import os
from typing import List, Tuple
from copy import deepcopy

from excitingtools.dict_utils import delete_nested_key

from ..exciting_settings.constants import keys_to_remove
from ..io.file_system import create_run_dir, copy_calculation_inputs, flatten_directory
from ..io.parsers import read_output_file
from ..io.tolerances import get_json_tolerances
from ..tester.compare import ErrorFinder
from ..tester.failure import Failure
from ..tester.report import TestResults, SummariseTests, skipped_test_summary

from .execute import execute


def remove_untested_keys(data: dict, full_file_name: str, keys_to_remove: dict) -> dict:
    """
    Remove keys not to be tested from data dictionary.

    :param dict data: Test or reference data, or Failure object
    :param str full_file_name: File name prepended by full path (relative to test/)
    :param dict keys_to_remove: dictionary of keys not to be tested
    :return dict data: Dictionary, with keys not meant for testing removed.
    """
    if not isinstance(data, dict):
        return data

    file_name = os.path.basename(full_file_name)

    # Remove all scf loops except for the last one
    if file_name == 'INFO.OUT':
        try:
            scl_indices = [int(item) for item in list(data['scl'].keys())]
            last_scl = max(scl_indices)
            data['scl'] = data['scl'][str(last_scl)]
        except KeyError:
            raise KeyError("Expected tolerances for INFO.OUT but no key 'scl' can be found.\n"
                           "Either tolerances are not for INFO.OUT, INFO.OUT's parser has changed, or "
                           "the 'scl' key has already been removed")

    # Remove untested keys
    try:
        for key_chain in keys_to_remove[file_name]:
            delete_nested_key(data, key_chain)
    except KeyError:
        return data

    return data


def compare_outputs(run_dir: str, ref_dir: str, output_files: List[str], tolerance: dict) -> dict:
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

        ref_data = remove_untested_keys(ref_data, file, keys_to_remove)
        test_data = remove_untested_keys(test_data, file, keys_to_remove)

        try:
            file_tolerances = tolerance[file]
        except KeyError:
            raise KeyError(f"File name is not a key in the tolerance file: {file}. \n"
                           "This may imply that tolerance template uses the wrong file name in its key.")

        file_results = ErrorFinder(test_data, ref_data, file_tolerances, file_name=file)
        test_results[file] = deepcopy(file_results)

    return test_results


def execute_and_compare_single_test(test_dir: str,
                                    full_run_dir: str,
                                    full_ref_dir: str,
                                    execute_cmd: str,
                                    main_out: str,
                                    max_time: int,
                                    output_files: List[str],
                                    just_tolerances: dict) -> TestResults:
    """
    Execute a test case and compare the output files to reference files.

    :param str test_dir: Path to test case (relative to test_farm).
    :param str full_run_dir: test_dir + reference_dir.
    :param str full_ref_dir: test_dir + run_dir.
    :param str execute_cmd: executable command.
    :param str main_out: File always output by program.
    :param int max_time: Time before program is terminated.
    :param List[str] output_files: List of files under test.
    :param dict just_tolerances: Tolerances without units.

    :return TestResults test_results: Test case results.
    """
    run_success, err_mess, timing = execute(full_run_dir, execute_cmd, main_out, max_time)

    test_results = TestResults(test_dir, run_success, timing)

    flatten_directory(full_run_dir)

    if run_success:
        test_results_dict = compare_outputs(full_run_dir, full_ref_dir, output_files, just_tolerances)
        test_results.set_results(test_results_dict)

    return test_results


def run_single_test(test_dir: str,
                    test_properties: dict,
                    run_dir: str,
                    ref_dir: str,
                    main_out: str,
                    executable: str,
                    max_time: int,
                    handle_errors: bool,
                    repeated_tests: dict) -> TestResults:
    """
    Runs a test case and compares output files under test against references.

    :param str test_dir:          Test name, prepended with relative path.
    :param dict test_properties:  Test properties, as defined by ConfigurationDefaults namedTuple.
    :param str main_out:          main output file of the exciting calculation
    :param str run_dir:           name of the run directory of the test case
    :param str ref_dir:           name of the ref directory of the test case
    :param str executable:        executable command
    :param int max_time:          max time a job is allowed to run for before being killed
    :param bool handle_errors:    Whether or not failures and skips are allowed to propagate
    :param dict repeated_tests:   (key:value) = (Test Name: N times to repeat)

    :return TestResults test_results: Test case results object.
    """
    method, test_name = test_dir.split('/')[-2:]
    full_ref_dir = os.path.join(test_dir, ref_dir)
    full_run_dir = os.path.join(test_dir, run_dir)

    print('Run test %s:' % test_name)

    tolerances_without_units, output_files = get_json_tolerances(full_ref_dir, test_properties['files_under_test'])
    create_run_dir(test_dir, run_dir)
    copy_calculation_inputs(full_ref_dir, full_run_dir, test_properties['inputs'])

    test_results = execute_and_compare_single_test(test_dir,
                                                   full_run_dir,
                                                   full_ref_dir,
                                                   executable,
                                                   main_out,
                                                   max_time,
                                                   output_files,
                                                   tolerances_without_units)
    test_results.print_results()

    # Rerun the test if it fails and is assigned as flakey in failing_tests.py
    if len(test_results.files_with_errors) > 0:
        try:
            n_repeats = repeated_tests[os.path.join(method, test_name)]
        except KeyError:
            n_repeats = 0

        for ith_repeat in range(1, n_repeats + 1):
            print(f'Repeating: {test_dir}: {ith_repeat} / {n_repeats}')
            test_results = execute_and_compare_single_test(test_dir,
                                                           full_run_dir,
                                                           full_ref_dir,
                                                           executable,
                                                           main_out,
                                                           max_time,
                                                           output_files,
                                                           tolerances_without_units)
            test_results.print_results()

            if len(test_results.files_with_errors) == 0:
                continue

    test_results.assert_errors(handle_errors)

    return test_results


def run_tests(main_out: str,
              tests: dict,
              run_dir: str,
              ref_dir: str,
              executable: str,
              np: int,
              omp: int,
              max_time: int,
              handle_errors: bool,
              repeated_tests: dict
              ) -> SummariseTests:
    """
    Runs tests in test_list.

    :param str main_out:          main output file of the exciting calculation
    :param dict tests:            All test cases that will be run
    :param str run_dir:           name of the run directory of the test case
    :param str ref_dir:           name of the ref directory of the test case
    :param str executable:        executable command for the exciting run
    :param int np:                number of MPI processes
    :param int omp:               number of OMP threads
    :param int max_time:          max time before a job is killed
    :param bool handle_errors:    Whether or not failures and skips are allowed to propagate
    :param dict repeated_tests:  (key:value) = (Test Name: N times to repeat)

    :return SummariseTests report: Test suite summary
    """
    if 'exciting_serial' in executable:
        print('Run tests with exciting_serial.')
    elif 'exciting_smp' in executable:
        print('Run tests with exciting_smp with %i open MP threads.' % (omp))
    elif 'exciting_mpismp' in executable:
        print('Run tests with exciting_mpismp with %i open MP threads and %i MPI processes.' % (omp, np))
    elif 'exciting_purempi' in executable:
        print('Run tests with exciting_purempi %i MPI processes.' % np)

    report = SummariseTests()
    for test_name, test_properties in tests.items():
        test_results = run_single_test(test_name,
                                       test_properties,
                                       run_dir,
                                       ref_dir,
                                       main_out,
                                       executable,
                                       max_time,
                                       handle_errors,
                                       repeated_tests)
        report.add(test_results)

    return report
