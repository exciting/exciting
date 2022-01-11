"""
Module containing:
* TestResults class, which handles the results of a single test case
* SummariseTests class, which summarises several TestResults instances
"""
from typing import List, Optional
import numpy as np

from .failure import Failure
from ..utilities.termcolor_wrapper import print_color
from ..tester.compare import ErrorFinder


class TestResults:
    def __init__(self, name: str, completed: bool, timing: float, test_results: Optional[dict] = None):
        """
        Initialise NewTest class to store the assertion results for a set of regression-tested files,
        for a given test case.

        :param str name: Test case name
        :param bool completed: Test case completed execution
        :param dict test_results: Dictionary of errors and tolerances for each tested file in the test case.
        """
        self.test_name = name
        self.completed = completed
        self.timing = timing

        self.file_names: List[str] = []
        self.n_files = 0
        self.results: List[ErrorFinder] = []
        self.succeeded_files: List[str] = []
        self.files_with_errors: List[str] = []
        self.files_with_io_error: List[str] = []
        # List of reference (or target) files that could not be opened/evaluated.
        self.unevaluated_files: List[str] = []

        if test_results is not None:
            value_types = [x for x in test_results.values()]
            assert isinstance(all(value_types), ErrorFinder), 'All entries in test_results should be class ErrorFinder'
            self.set_results(test_results)

    def set_results(self, test_results: dict):
        if self.completed:
            # All results
            self.results = [file_result for file_result in test_results.values()]
            self.n_files = len(self.results)
            # Unevaluated files
            self.unevaluated_files = [name for name, file_result in test_results.items()
                                      if isinstance(file_result, Failure)]
            # Evaluated files
            evaluated_results = [file_result for name, file_result in test_results.items()
                                 if not isinstance(file_result, Failure)]

            self.file_names = np.asarray([file_result.file_name for file_result in evaluated_results])
            n_errors_per_file = np.asarray([file_result.n_errors for file_result in evaluated_results])
            # List of files names for which all comparisons succeeded
            self.succeeded_files = self.file_names[np.where(n_errors_per_file == 0)]
            # List of files for which any comparison failed
            self.files_with_errors = self.file_names[np.where(n_errors_per_file != 0)]

    def print_results(self):
        """
        Wrapper function that prints a header for the test case, then iterates over
        each file comparison and calls its print method.

        For example, print the header string:
           Test Case: test_farm/properties/LDA_PZ-wannier-SiC/run
           Time (s): 13.6
           Report: SUCCESS 8/9, UNEVALUATED 0/9, FAIL 1/9.

        Followed by the result of each compared file
        """
        print('Test Case:', self.test_name)
        print('Time (s): %.1f' % self.timing)

        if not self.completed:
            print('Test execution failed to complete')
            return

        assert self.results is not None, "results in TestResults must be filled before calling print_results"

        print('Test execution completed')
        print('Output Files:', end=' ')
        print_color('SUCCESS %i/%i' % (len(self.succeeded_files), self.n_files), 'green', end='')
        print(',', end=' ')
        print_color('NOT EVALUATED %i/%i' % (len(self.unevaluated_files), self.n_files), 'yellow', end='')
        print(',', end=' ')
        print_color('FAIL %i/%i' % (len(self.files_with_errors), self.n_files), 'red', end='')
        print('.')

        for file_result in self.results:
            file_result.print_result()

    def assert_errors(self, handle_errors: bool, fail_if_unevaluated=True):
        """
        Throws if the any assertions for if a test case failed.

        Failure is defined as any comparison/assertion failing for a given file, or if a reference/target file
        could not be opened (implying either the exciting run did not write it, or the reference is missing).

        :param bool handle_errors: If true, allow any errors to propagate to the end of the test suite
        :param bool fail_if_unevaluated: If true, unevaluated files count as a failure
        """
        if not handle_errors:
            if fail_if_unevaluated:
                assert (len(self.files_with_errors) == 0) and (len(self.unevaluated_files) == 0), "Test case failed"
            else:
                assert len(self.files_with_errors) == 0, "Test case failed"


def timing_statistics(timing_dict: dict) -> dict:
    """
    Test suite timings

    :param dict timing_dict: Timing dictionary of the form {'name': time}
    :return dict timing_results: Results for timings
    """

    # Create a list of timings and a list of names, with same ordering
    names = []
    times = []
    for name, time in timing_dict.items():
        names.append(name)
        times.append(time)

    total_time = sum(times)
    avg_time = (total_time / len(times)) if len(times) > 0 else 0

    timing_results = {'total': total_time,
                      'average': avg_time,
                      'shortest': {'name': 'null', 'time': 0},
                      'longest': {'name': 'null', 'time': 0}
                      }

    # Null run
    if len(times) == 0:
        return timing_results

    # Shortest time
    i_min = np.argmin(times)
    timing_results['shortest'] = {'name': names[i_min], 'time': times[i_min]}

    # Longest time
    i_max = np.argmax(times)
    timing_results['longest'] = {'name': names[i_max], 'time': times[i_max]}

    return timing_results


class SummariseTests:
    def __init__(self):
        """
        Summarise the test suite results.
        Individual reporting is done by class TestResults
        """
        # Test cases
        self.n_test_cases = 0
        self.n_succeeded_test_cases = 0
        self.n_failed_test_cases = 0
        self.n_unevaluated_files = 0

        # Files (over all test cases)
        self.n_files = 0
        self.n_succeeded_files = 0
        self.n_failed_files = 0

        # Timing
        self.times = {}

    def add(self, test: TestResults):
        """
        Add results from a test case (comprised of regression-testing on multiple output files)
        """
        # Test Case
        errors_per_file = np.asarray([file_result.n_errors for file_result in test.results
                                      if not isinstance(file_result, Failure)])
        successful_files = np.where(errors_per_file == 0)[0]
        failed_files = np.where(errors_per_file != 0)[0]

        self.n_files += len(test.results)
        self.n_succeeded_files += len(successful_files)
        self.n_failed_files += len(failed_files)
        self.n_unevaluated_files += len(test.unevaluated_files)
        assert (self.n_succeeded_files + self.n_failed_files + self.n_unevaluated_files) == self.n_files

        # Test Suite
        succeeded_test = len(successful_files) == len(test.results)
        self.n_test_cases += 1
        self.n_succeeded_test_cases += int(succeeded_test)
        self.n_failed_test_cases += int(not succeeded_test)

        self.times.update({test.test_name: test.timing})

    def print(self):
        """
        Print summary of multiple test cases.
        """
        summary_msg = f"Test Suite Summary\n\n"

        summary_msg += f"Total number of files: {self.n_files}\n"
        summary_msg += f"Succeeded file comparisons: {self.n_succeeded_files}\n"
        summary_msg += f"Failed file comparisons: {self.n_failed_files}\n"
        summary_msg += f"Unevaluated files: {self.n_unevaluated_files}\n\n"

        summary_msg += f"Total number of test cases: {self.n_test_cases}\n"
        summary_msg += f"Succeeded test cases: {self.n_succeeded_test_cases}\n"
        summary_msg += f"Failed test cases: {self.n_failed_test_cases}\n"
        # TODO(Alex). Issue 98. Move Printing of Skipped Test Cases
        # Currently handled by a free function, using the failing_tests.py file.
        # summary_msg += f"Skipped test cases: ADD ME \n\n "
        print(summary_msg)

    def print_timing(self):
        """
        Print timing summary for the test suite
        """
        timing_results = timing_statistics(self.times)
        total_time = timing_results['total']

        unit_str = 's'
        if total_time >= 60:
            total_time = total_time / 60.
            unit_str = 'mins'

        print('Total test suite time (' + unit_str + ') : %.1f' % total_time)
        print('Average test time (s): %.1f' % timing_results['average'])
        print('Shortest test time (s): %.1f, taken by %s' %
              (timing_results['shortest']['time'], timing_results['shortest']['name']))
        print('Longest test time (s): %.1f, taken by %s' %
              (timing_results['longest']['time'], timing_results['longest']['name']))


def skipped_test_summary(skipped_tests: dict, print_comments=False):
    """
    Summarises skipped tests

    :param list skipped_tests: list of skipped tests.
           Each element is a dictionary of the skipped test properties
    :param bool print_comments: Whether or not to print failing comments associated with failing tests.
    :return None
    """
    if len(skipped_tests) == 0:
        return

    if print_comments:
        print_func = lambda name, properties: print(' ', name, '. ', properties['comments'])
    else:
        print_func = lambda name, properties: print(' ', name, '. ')

    print('Summary of SKIPPED tests (as defined in tests_config.yml):')
    for name, properties in skipped_tests.items():
        print_func(name, properties)


def indent(elem, level=0):
    i = "\n" + level * "  "
    j = "\n" + (level - 1) * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for subelem in elem:
            indent(subelem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = j
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = j
    return elem
