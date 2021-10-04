"""
Test the class SummariseTests attributes and print method

Running
-----------
Should be run from exciting/test. From exciting's root:
```
cd test
pytest -s test_functions/test_class_summarise_tests.py
```
else the relative imports will cause test failures

Debugging Notes
-----------
Debugging differences in string buffers can be difficult due to inconsistent whitespaces. If one runs with:
  `pytest -svv test_functions/test_run_single_test.py`
pytest will output exactly what is compared in the assert, so one can copy the strings and check in
(for example) a text editor that they're entirely consistent.

Note, commenting out and subsequently uncommenting the expected buffer strings *can* modify the formatting
w.r.t. whitespaces. Maybe this method of testing is not robust.

TODO
Test timings
Test unevaluated_files
"""
import pytest
import sys
from io import StringIO
import os
from typing import List

# Rename TestResults else pytest will try to collect from it
from ..modules.tester.report import TestResults as ExTestResults, SummariseTests
from ..modules.runner.test import load_tolerances, strip_tolerance_units, compare_outputs_json


def redirect_stdout(func):
    """
    Decorator to redirect stdout to a string.
    :param func: A print function
    """

    def wrapper():
        old_stdout = sys.stdout
        sys.stdout = my_stdout = StringIO()
        func()
        sys.stdout = old_stdout
        return my_stdout.getvalue()

    return wrapper


@pytest.fixture()
def successful_test_dir() -> str:
    return "test_functions/dummy_app_tests/LDA_VWN-He_passing"


@pytest.fixture()
def failing_arrays1d2d_test_dir() -> str:
    return "test_functions/dummy_app_tests/GW_ZrO2_failing_arrays"


def get_test_results(test_dir: str, output_files_to_check: List[str]) -> dict:
    """
    Set the run and reference directories, files to regression test, tolerances, and
    return the results of comparing the test data.
    """
    full_run_dir = os.path.join(test_dir, "run")
    full_ref_dir = os.path.join(test_dir, "ref")

    json_tolerances = strip_tolerance_units(load_tolerances(full_ref_dir))
    test_results_dict = compare_outputs_json(full_run_dir, full_ref_dir, output_files_to_check, json_tolerances)
    return test_results_dict


def test_print_successful(successful_test_dir):

    gs_files = ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml']
    test_results = ExTestResults(successful_test_dir, completed=True, timing=0.0)
    test_results_dict = get_test_results(successful_test_dir, gs_files)
    test_results.set_results(test_results_dict)

    report = SummariseTests()
    report.add(test_results)

    assert report.n_test_cases == 1, "Expect 1 test case"
    assert report.n_files == 5, "Comprised of 5 file comparisons (see gs_files)"

    assert report.n_succeeded_test_cases == 1, "Test case passes"
    assert report.n_failed_test_cases == 0,  "No test cases fail"
    assert report.n_succeeded_files == 5, "All files should pass"
    assert report.n_failed_files == 0, "No files should fail"
    assert report.n_unevaluated_files == 0, "All files should be evaluated"

    expected_print_buffer = """Test Suite Summary

Total number of files: 5
Succeeded file comparisons: 5
Failed file comparisons: 0
Unevaluated files: 0

Total number of test cases: 1
Succeeded test cases: 1
Failed test cases: 0

"""

    @redirect_stdout
    def print_results():
        return report.print()

    output_print_buffer = print_results()

    assert output_print_buffer == expected_print_buffer, \
        "Expect the print() method of SummariseTests to give expected_print_buffer"


def test_print_failing_testcase(failing_arrays1d2d_test_dir):

    test_results = ExTestResults(failing_arrays1d2d_test_dir, completed=True, timing=0.0)
    test_results_dict = get_test_results(failing_arrays1d2d_test_dir, ['EVALQP.DAT'])
    test_results.set_results(test_results_dict)

    report = SummariseTests()
    report.add(test_results)

    assert report.n_test_cases == 1, "Expect 1 test case"
    assert report.n_files == 1, "Comprised of 1 file comparison (EVALQP.DAT)"

    assert report.n_succeeded_test_cases == 0, "No test cases pass"
    assert report.n_failed_test_cases == 1,  "One test case fails"
    assert report.n_succeeded_files == 0, "No files should pass"
    assert report.n_failed_files == 1, "The only file should fail"
    assert report.n_unevaluated_files == 0, "All files should be evaluated"

    @redirect_stdout
    def print_results():
        return report.print()

    output_print_buffer = print_results()

    expected_print_buffer = """Test Suite Summary

Total number of files: 1
Succeeded file comparisons: 0
Failed file comparisons: 1
Unevaluated files: 0

Total number of test cases: 1
Succeeded test cases: 0
Failed test cases: 1

"""

    assert output_print_buffer == expected_print_buffer, \
        "Expect the print() method of SummariseTests to give expected_print_buffer"


def test_timings(successful_test_dir):
    """
    Test the SummariseTests class for timings
    """
    # Initialised with no test results, other than name, completed and timing
    slow_test = ExTestResults(name='Slow test', completed=True, timing=600.0)
    avg_test = ExTestResults(name='Average-speed test', completed=True, timing=10.0)
    fast_test = ExTestResults(name='Fast test', completed=True, timing=0.99)

    report = SummariseTests()
    report.add(slow_test)
    report.add(avg_test)
    report.add(fast_test)

    @redirect_stdout
    def print_times():
        return report.print_timing()

    output_print_buffer = print_times()

    expected_buffer = """Total test suite time (mins) : 10.2
Average test time (s): 203.7
Shortest test time (s): 1.0, taken by Fast test
Longest test time (s): 600.0, taken by Slow test
"""

    assert output_print_buffer == expected_buffer, \
        "Expect the print_timing() method of SummariseTests to give expected_buffer"
