"""
Test the classes TestResults

Running
-----------
Should be runnable from any folder
```
cd test
pytest -s test_functions/test_class_testresults.py
```

Debugging Notes
-----------
Debugging differences in string buffers can be difficult due to inconsistent whitespaces. If one runs with:
  `pytest -svv test_functions/test_run_single_test.py`
pytest will output exactly what is compared in the assert, so one can copy the strings and check in
(for example) a text editor that they're entirely consistent.

Note, commenting out and subsequently uncommenting the expected buffer strings *can* modify the formatting
w.r.t. whitespaces. Maybe this method of testing is not robust.

TODO(A/B/H) Issue 96. Test Error Logging with 3D Array
"""
import pytest
import sys
from io import StringIO
import os
from typing import List, Union
import re
from pathlib import Path

from ..src.tester.compare import ErrorFinder
# Rename TestResults else pytest will try to collect from it
from ..src.tester.report import TestResults as ExTestResults
from ..src.io.tolerances import load_tolerances, strip_tolerance_units, list_tolerance_files_in_directory
from ..src.runner.test import compare_outputs


def strip_ansi(text: str) -> str:
    r"""
    Removes all ansi directives from the string.

    References:
        http://stackoverflow.com/questions/14693701/remove-ansi
        https://stackoverflow.com/questions/13506033/filtering-out-ansi-escape-sequences

    Examples:
        line = '\t\u001b[0;35mBlabla\u001b[0m     \u001b[0;36m172.18.0.2\u001b[0m'
        escaped_line = strip_ansi(line)
        assert escaped_line == '\tBlabla     172.18.0.2'
    """
    ansi_escape3 = re.compile(r'(\x9B|\x1B\[)[0-?]*[ -/]*[@-~]', flags=re.IGNORECASE)
    return ansi_escape3.sub('', text)


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


def get_test_framework_root():
    this_file_path = os.path.dirname(os.path.realpath(__file__))
    exciting_root = Path(this_file_path).parents[1]
    return exciting_root / 'test'


@pytest.fixture()
def successful_test_dir() -> str:
    return "test_functions/dummy_app_tests/LDA_VWN-He_passing"


@pytest.fixture()
def failing_scalars_test_dir() -> str:
    return "test_functions/dummy_app_tests/LDA_VWN-He_failing_scalars"


@pytest.fixture()
def failing_arrays1d2d_test_dir() -> str:
    return "test_functions/dummy_app_tests/GW_ZrO2_failing_arrays"


@pytest.fixture()
def missing_files_test_dir() -> Path:
    root = get_test_framework_root()
    return root / "test_functions/dummy_app_tests/LDA_VWN-He_missing_files"


def get_test_results(test_dir: Union[str, Path], output_files_to_check: List[str]) -> dict:
    """
    Set the run and reference directories, files to regression test, tolerances, and
    return the results of comparing the test data.
    """
    root = get_test_framework_root()
    full_run_dir = os.path.join(root, test_dir, "run")
    full_ref_dir = os.path.join(root, test_dir, "ref")

    tolerance_files = list_tolerance_files_in_directory(full_ref_dir)
    json_tolerances = load_tolerances(full_ref_dir, tolerance_files)
    json_tolerances = strip_tolerance_units(json_tolerances)

    test_results_dict = compare_outputs(full_run_dir, full_ref_dir, output_files_to_check, json_tolerances)
    return test_results_dict


def test_attributes_successful(successful_test_dir):
    """
    Test the attributes of class TestResults, for a successful test case.
    """
    gs_files = ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml']
    test_results = ExTestResults(successful_test_dir, completed=True, timing=0.0)
    test_results_dict = get_test_results(successful_test_dir, gs_files)
    test_results.set_results(test_results_dict)

    assert os.path.basename(test_results.test_name) == "LDA_VWN-He_passing", "Base test name"
    assert test_results.completed, "exciting run completed (rather than crashed)"
    assert set(test_results.file_names) == set(gs_files), "Expect these ground state files to be regression-tested"
    assert set(test_results.succeeded_files) == set(gs_files), "All tested files should succeed"
    assert len(test_results.files_with_errors) == 0, "All tested files should succeed"
    assert len(test_results.unevaluated_files) == 0, "All listed files should be tested"


def test_print_successful(successful_test_dir):
    """
    Test TestResults print() method on a successful test case.
    """
    gs_files = ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml']
    test_results = ExTestResults(successful_test_dir, completed=True, timing=0.0)
    test_results_dict = get_test_results(successful_test_dir, gs_files)

    for name, file_results in test_results_dict.items():
        message = "If not ErrorFinder, file_results returns Failure => " + name + " missing"
        assert isinstance(file_results, ErrorFinder), message

    test_results.set_results(test_results_dict)

    @redirect_stdout
    def print_results():
        return test_results.print_results()

    output_print_buffer = print_results()
    colourless_output_buffer = strip_ansi(output_print_buffer)

    expected_print_buffer = f"""Test Case: {test_results.test_name}
Time (s): 0.0
Test execution completed
Output Files: SUCCESS 5/5, NOT EVALUATED 0/5, FAIL 0/5.
    INFO.OUT SUCCESS
    evalcore.xml SUCCESS
    geometry.xml SUCCESS
    eigval.xml SUCCESS
    atoms.xml SUCCESS
"""
    assert colourless_output_buffer == expected_print_buffer, \
        "Expect the print() method of TestResults to give expectd_print_buffer"


def test_attributes_erroneous_scalars(failing_scalars_test_dir):
    """
    Test the TestResults class attributes when some of the run data is defined to be slightly different
    to the reference data.
    """
    gs_files = ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml']
    test_results_dict = get_test_results(failing_scalars_test_dir, gs_files)

    for name, file_results in test_results_dict.items():
        message = "If not ErrorFinder, file_results returns Failure => " + name + " missing"
        assert isinstance(file_results, ErrorFinder), message

    test_results = ExTestResults(failing_scalars_test_dir, completed=True, timing=0.0)
    test_results.set_results(test_results_dict)

    assert os.path.basename(test_results.test_name) == "LDA_VWN-He_failing_scalars", "Base test name"
    assert test_results.completed, "exciting run completed (rather than crashed)"
    assert set(test_results.file_names) == set(gs_files), "Expect specified ground state files to be regression-tested"
    assert set(test_results.succeeded_files) == {'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml'}, \
        "All files should succeed apart from INFO.OUT"
    assert test_results.files_with_errors == ["INFO.OUT"], "run/INFO.OUT should contain discrepancies w.r.t. ref/INFO.OUT"
    assert len(test_results.unevaluated_files) == 0, "All listed files in `gs_files` should be tested"


def test_print_erroneous_scalars(failing_scalars_test_dir):
    """
    Test TestResults print() method on a failing test case
    """
    gs_files = ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml']
    test_results_dict = get_test_results(failing_scalars_test_dir, gs_files)

    for name, file_results in test_results_dict.items():
        message = "If not ErrorFinder, file_results returns Failure => " + name + " missing"
        assert isinstance(file_results, ErrorFinder), message

    test_results = ExTestResults(failing_scalars_test_dir, completed=True, timing=0.0)
    test_results.set_results(test_results_dict)

    @redirect_stdout
    def print_results():
        return test_results.print_results()

    output_print_buffer = print_results()
    colourless_output_buffer = strip_ansi(output_print_buffer)
    
    expected_print = f"""Test Case: {test_results.test_name}
Time (s): 0.0
Test execution completed
Output Files: SUCCESS 4/5, NOT EVALUATED 0/5, FAIL 1/5.
    INFO.OUT FAIL
      Failures TOTAL: 3, INTEGER: 1, FLOAT: 1, ARRAY: 0, STRING 1
      INTEGER FAILURES
      Key                                          Test Data      Ref Data          Diff     Tolerance
      --------------------------------------------------------------------------------------------------------------
      initialization%APW functions                         8             4             4             0

      FLOAT FAILURES
      Key                                          Test Data      Ref Data          Diff     Tolerance
      --------------------------------------------------------------------------------------------------------------
      scl%Total energy                           -2.83483618  -20.83483618      1.80e+01      1.00e-06

      STRING FAILURES
      Key                                                               Test Data                           Ref Data
      --------------------------------------------------------------------------------------------------------------
      initialization%Spin treatment                              spin-unpolarised                     spin-polarised

    evalcore.xml SUCCESS
    geometry.xml SUCCESS
    eigval.xml SUCCESS
    atoms.xml SUCCESS
"""

    assert colourless_output_buffer == expected_print, \
        "Expect the print() method of TestResults to give expected_print_buffer"


def test_attributes_failing_array(failing_arrays1d2d_test_dir):
    """
    Test the TestResults class attributes when some of the run data is defined to be different beyond tolerance
    to the reference data.
    """
    gw_files = ['EVALQP.DAT']
    test_results_dict = get_test_results(failing_arrays1d2d_test_dir, gw_files)

    for name, file_results in test_results_dict.items():
        message = "If not ErrorFinder, file_results returns Failure => " + name + " missing"
        assert isinstance(file_results, ErrorFinder), message

    test_results = ExTestResults(failing_arrays1d2d_test_dir, completed=True, timing=0.0)
    test_results.set_results(test_results_dict)

    assert os.path.basename(test_results.test_name) == "GW_ZrO2_failing_arrays", "Base test name"
    assert test_results.completed, "exciting run completed (rather than crashed)"
    assert set(test_results.file_names) == {'EVALQP.DAT'}, "Expect only this GW file to be regression-tested"
    assert len(test_results.succeeded_files) == 0, "Quasi-particle run file should contain errors"
    assert test_results.files_with_errors == ["EVALQP.DAT"], "run/EVALQP.DAT should contain discrepancies w.r.t. ref/EVALQP.DAT"
    assert len(test_results.unevaluated_files) == 0, "All listed files should be tested"


def test_print_failing_array(failing_arrays1d2d_test_dir):
    """
    Test TestResults print() method on a successful test case.
    """
    gw_files = ['EVALQP.DAT']
    test_results_dict = get_test_results(failing_arrays1d2d_test_dir, gw_files)

    for name, file_results in test_results_dict.items():
        message = "If not ErrorFinder, file_results returns Failure => " + name + " missing"
        assert isinstance(file_results, ErrorFinder), message

    test_results = ExTestResults(failing_arrays1d2d_test_dir, completed=True, timing=0.0)
    test_results.set_results(test_results_dict)

    @redirect_stdout
    def print_results():
        return test_results.print_results()

    output_print_buffer = print_results()
    colourless_output_buffer = strip_ansi(output_print_buffer)

    expected_print_buffer = f"""Test Case: {test_results.test_name}
Time (s): 0.0
Test execution completed
Output Files: SUCCESS 0/1, NOT EVALUATED 0/1, FAIL 1/1.
    EVALQP.DAT FAIL
      Failures TOTAL: 3, INTEGER: 0, FLOAT: 0, ARRAY: 3, STRING 0
      ARRAY FAILURES
      Key: 1%k_point
      (element)                                    Test Data      Ref Data          Diff     Tolerance
      --------------------------------------------------------------------------------------------------------------
      (0)                                         0.10000000    0.00000000      1.00e-01      1.00e-08

      Key: 1%energies
      (element)                                    Test Data      Ref Data          Diff     Tolerance
      --------------------------------------------------------------------------------------------------------------
      (9,3)                                      -1.52291000   -1.52299000      8.00e-05      1.00e-06
      (23,6)                                     -0.47984000   -0.47985000      1.00e-05      1.00e-06

      Key: 2%energies
      (element)                                    Test Data      Ref Data          Diff     Tolerance
      --------------------------------------------------------------------------------------------------------------
      (18,2)                                     -0.13826000   -0.19826000      6.00e-02      1.00e-06
      (49,7)                                      0.34371000    0.34271000      1.00e-03      1.00e-06

"""

    assert colourless_output_buffer == expected_print_buffer, \
        "Expect the print() method of TestResults to give expected_print_buffer"


def test_attributes_missing_files(missing_files_test_dir):
    """
    Test the TestResults class attributes when a file is missing
    """
    gs_files = ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml']
    present_files = ['evalcore.xml', 'eigval.xml']
    missing_files = ['INFO.OUT', 'atoms.xml', 'geometry.xml']
    test_results_dict = get_test_results(missing_files_test_dir, gs_files)

    test_results = ExTestResults(missing_files_test_dir, completed=True, timing=0.0)
    test_results.set_results(test_results_dict)

    assert os.path.basename(test_results.test_name) == "LDA_VWN-He_missing_files", "Base test name"
    assert test_results.completed, "exciting run completed (rather than crashed)"
    assert set(test_results.file_names) == set(present_files), "Expect present files to be regression-tested"
    assert set(test_results.succeeded_files) == set(present_files), "All files that do not give an IO error should succeed"
    assert len(test_results.files_with_errors) == 0, "All tested files should succeed"
    assert set(test_results.unevaluated_files) == set(missing_files), "All listed files except for INFO.OUT and atoms.xml should be tested"


def test_print_missing_file(missing_files_test_dir):
    """
    Test TestResults print() method when a file is missing
    """
    gs_files = ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml']
    test_results = ExTestResults(missing_files_test_dir, completed=True, timing=0.0)
    test_results_dict = get_test_results(missing_files_test_dir, gs_files)

    test_results.set_results(test_results_dict)

    @redirect_stdout
    def print_results():
        return test_results.print_results()

    output_print_buffer = print_results()
    colourless_output_buffer = strip_ansi(output_print_buffer)
    root = get_test_framework_root()

    expected_print_buffer = f"""Test Case: {test_results.test_name}
Time (s): 0.0
Test execution completed
Output Files: SUCCESS 2/5, NOT EVALUATED 3/5, FAIL 0/5.
    INFO.OUT NOT EVALUATED
      File {root}/test_functions/dummy_app_tests/LDA_VWN-He_missing_files/run/INFO.OUT could not be opened.
    evalcore.xml SUCCESS
    geometry.xml NOT EVALUATED
      File {root}/test_functions/dummy_app_tests/LDA_VWN-He_missing_files/ref/geometry.xml could not be opened.
    eigval.xml SUCCESS
    atoms.xml NOT EVALUATED
      File {root}/test_functions/dummy_app_tests/LDA_VWN-He_missing_files/run/atoms.xml could not be opened.
"""

    assert colourless_output_buffer == expected_print_buffer, \
        "Expect the print() method of TestResults to give expected_print_buffer"
