import os
from typing import List, Union, Tuple
import json
from copy import deepcopy
import re

from excitingtools.dict_utils import delete_nested_key

from ..constants import settings, methods_moved_to_json, keys_to_remove
from ..parsers import ErroneousFileError, parser_chooser
from ..infrastructure import create_run_dir, copy_exciting_input, get_test_from_init, flatten_directory
from ..wildcard_processor import wildcard_processor

from ..termcolor_wrapper import print_color
from ..tester.compare import ErrorFinder
from ..tester.test import from_init
from ..tester.report import Report, test_suite_summary, skipped_test_summary, timing_summary, TestResults, SummariseTests
from ..tester.failure import Failure, Failure_code

from .execute import execute
from ..wildcard_processor import wildcard


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
        return data
    except OSError:
        failure_code = {'ref': Failure_code.REFERENCE, 'run': Failure_code.RUN}
        return Failure(test_name=file_name, failure_code=failure_code[sub_dir])
    except ErroneousFileError:
        return Failure(test_name=file_name, failure_code=Failure_code.ERRORFILE)


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


# TODO(Alex) Ultimately we can remove this and have the config file deal with it (in a few merges)
def update_wildcard_files_under_test(files_under_test: List[str], file_names: List[str]) -> List[str]:
    """
    Create copies of 'files_under_test' entries for every file name with a common name prefix.

    For files for which there can be multiple extensions (FILE_0.OUT, FILE_1.OUT, etc), it is expected that
    a single tolerance will be defined with the key 'FILE_<WILDCARD>' in the tolerances dictionary and therefore
    only a single entry 'FILE_<WILDCARD>' in the 'files_under_test' list

    This routine will copy this element for every reference file that matches the prefix.

    For example:
      On input, files_under_test = ['FILE_<WILDCARD>', 'OTHER_XY.OUT']
      Reference files = ['FILE_ABC_0.OUT', 'FILE_DEF_1.OUT', 'OTHER_XY.OUT']
      On output, files_under_test_updated = ['FILE_ABC_0.OUT', 'FILE_DEF_1.OUT', 'OTHER_XY.OUT']

    :param List[str] files_under_test: File names under regression testing, which may or may not contain wildcards in their
     names.
    :param List[str] file_names: List of file names containing (but not exclusively containing) all reference output
    files to test i.e. ['INFO.OUT', 'GW_INFO.OUT', 'random_file.out'].
    :return List[str] files_under_test_updated: File names under regression testing, with wildcards replaced according
    to reference files listed in file_names.
    """
    assert files_under_test, "No 'files_under_test' are specified"
    files_under_test_updated = []
    for file_under_test in files_under_test:
        regex = re.compile(wildcard_processor(file_under_test))

        for file_name in file_names:
            if regex.match(file_name):
                files_under_test_updated.append(file_name)

    return files_under_test_updated


def update_wildcard_tolerances(tolerances: dict, files_under_test: List[str]) -> dict:
    """
    Create copies of tolerance entries for every file with a common name prefix.

    For files for which there can be multiple extensions (FILE_0.OUT, FILE_1.OUT, etc), it is expected that
    a single tolerance will be defined with the key 'FILE_<WILDCARD>' in the tolerances dictionary.

    This routine will copy this key:value for every reference file that matches the prefix.

    For example:
      On input, tolerances = {'FILE_<WILDCARD>': tol_dict}
      Reference file_names = [FILE_0.OUT, FILE_1.OUT]
      On output, tolerances = {'FILE_0.OUT': tol_dict, 'FILE_1.OUT': tol_dict,}

    If a key is already present, the initial tolerances will be retained.

    :param dict tolerances: Tolerances for all files under regression testing.
    :param List[str] files_under_test: List of all file names under testing (should not contain wildcards)
    :return dict tolerances_updated: Tolerances with for all files under regression testing.
    """
    for file in files_under_test:
        if any([True for w in wildcard if w in file]):
            raise ValueError('Wild card string present in file name: ' + file)

    tolerances_updated = {}
    unmatched_file_names = []

    for file_name, file_tolerances in tolerances.items():
        regex = re.compile(wildcard_processor(file_name))
        matched = False

        for reference_file in files_under_test:
            if regex.match(reference_file):
                tolerances_updated[reference_file] = file_tolerances
                matched = True

        if not matched:
            unmatched_file_names.append(file_name)

    if unmatched_file_names:
        raise KeyError(f"One or more keys in the tolerances JSON does not correspond to reference file names or "
                       f"cannot be matched with the wildcard symbols defined in wildcards.py: "
                       f"{unmatched_file_names}")

    return tolerances_updated


def list_tolerance_files_in_directory(directory: str) -> List[str]:
    """
    List the tolerances from any files in 'directory' of the form `*tol*.json`

    :param str directory: Directory containing files
    :return List[str] tolerance_files: List of tolerance files
    """
    files_in_directory = next(os.walk(directory))[2]
    r = re.compile("tolerance*.*json")
    tolerance_files = list(filter(r.match, files_in_directory))
    return tolerance_files


def load_tolerances(directory: str, tolerance_files: List[str]) -> dict:
    """
    Load JSON tolerance files

    :param str directory: Directory containing `*tol*.json` file
    :param List[str] tolerance_files:
    :return dict tolerances: Dictionary of tolerances
    """
    if not tolerance_files:
        raise FileNotFoundError('No tolerance files listed in `tolerance_files` argument')

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


def get_json_tolerances(full_ref_dir: str) -> Tuple[dict, List[str]]:
    """
    Load and preprocess JSON tolerance files.

    TODO(A/B/H) Issue 100. Update ErrorFinder class to use tol AND units in data comparison
       This would be better than to stripping the units (could pass the comparison function as an argument)
       See function documentation for a description of the issue

    :param str full_ref_dir: Full path to reference directory
    :return Tuple[dict, List[str]]: tolerances_without_units, reference_outputs: Tolerances dictionary
     of the form {'method': tols}, and output files to compare.
    """
    tolerance_files: List[str] = list_tolerance_files_in_directory(full_ref_dir)
    tolerances: dict = load_tolerances(full_ref_dir, tolerance_files)

    files_in_directory = next(os.walk(full_ref_dir))[2]
    reference_outputs = update_wildcard_files_under_test(tolerances.pop('files_under_test'), files_in_directory)

    tolerances = update_wildcard_tolerances(tolerances, reference_outputs)
    tolerances_without_units = strip_tolerance_units(tolerances)

    return tolerances_without_units, reference_outputs


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

        ref_data = remove_untested_keys(ref_data, file, keys_to_remove)
        test_data = remove_untested_keys(test_data, file, keys_to_remove)

        try:
            file_tolerances = tolerance[file]
        except KeyError:
            raise KeyError("File name does not match the file_name_key in tolerance file. "
                           "Most likely the key defined in tolerance templates is not equal "
                           "to the output file name")

        file_results = ErrorFinder(test_data, ref_data, file_tolerances, file_name=file)
        test_results[file] = deepcopy(file_results)

    return test_results


# TODO(Alex) Issue 99. Function to remove once test suite migrates to JSON
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
                         species_files: List[str], executable: str, max_time: int, handle_errors: bool,
                         repeated_tests: dict) -> TestResults:
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
    :param dict repeated_tests: (key:value) = (Test Name: N times to repeat)

    :return TestResults test_results: Test case results object.
    """
    head_tail = os.path.split(test_dir)
    test_name = head_tail[1]
    method = head_tail[0].split('/')[-1]
    full_ref_dir = os.path.join(test_dir, ref_dir)
    full_run_dir = os.path.join(test_dir, run_dir)

    print('Run test %s:' % test_name)

    tolerances_without_units, output_files = get_json_tolerances(full_ref_dir)
    create_run_dir(test_dir, run_dir)
    copy_exciting_input(full_ref_dir, full_run_dir, species_files, input_file)

    test_results = execute_and_compare_single_test_json(test_dir,
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
            test_results = execute_and_compare_single_test_json(test_dir,
                                                                full_run_dir,
                                                                full_ref_dir,
                                                                executable,
                                                                main_out,
                                                                max_time,
                                                                output_files,
                                                                just_tolerances)
            test_results.print_results()

            if len(test_results.files_with_errors) == 0:
                continue

    test_results.assert_errors(handle_errors)

    return test_results


def execute_and_compare_single_test_json(test_dir: str,
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
        test_results_dict = compare_outputs_json(full_run_dir, full_ref_dir, output_files, just_tolerances)
        test_results.set_results(test_results_dict)

    #test_results.print_results()
    #test_results.assert_errors(handle_errors)

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


def run_tests(main_out: str,
              test_list: List[str],
              run_dir: str,
              ref_dir: str,
              input_file: str,
              species_files: List[str],
              init_default: str,
              executable: str,
              np: int, omp: int,
              max_time: int,
              skipped_tests: List[str],
              handle_errors: bool,
              repeated_tests: dict
              ):
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
    :param dict repeated_tests: (key:value) = (Test Name: N times to repeat)
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
                                            executable, max_time, handle_errors, repeated_tests)
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
