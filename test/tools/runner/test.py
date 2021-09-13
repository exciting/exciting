import os
from typing import List

from ..parsers import ErroneousFileError, parser_chooser
from ..infrastructure import create_run_dir, copy_exciting_input, get_test_from_init, flatten_directory

from ..tester.test import from_init
from ..tester.report import Report, test_suite_summary, skipped_test_summary, timing_summary
from ..tester.failure import Failure, Failure_code

from .execute import execute


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
