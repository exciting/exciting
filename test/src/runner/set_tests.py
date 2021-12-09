"""
Process full test list to remove those that should not be run
for a given configuration
"""
import os
import re
from typing import List, Tuple

from ..exciting_settings.constants import settings
from ..exciting_settings.failing_tests import hanging_tests, failing_tests, repeat_tests
from ..io.parsers import get_compiler_type
from ..runner.profile import Compiler, Build_type, build_type_str_to_enum


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


def set_skipped_tests(executable: str, incl_failing_tests: bool) -> list:
    """
    Generate full list of tests to skip from the test suite

    :param str executable:
    :param bool incl_failing_tests: bool. If true, include failing tests in the suite,
    hence exclude from the list of skipped tests

    :return list tests_to_skip: list of tests to skip, where each entry is a dict
    """

    compiler = get_compiler_type()

    if compiler is None:
        print('Compiler could not be determined from build/make.inc\n'
              'Excluding all failing tests from the test suite')
        compiler = Compiler.all

    build_type = build_type_str_to_enum[executable]
    tests_to_skip = []

    # Always include hanging tests in tests to skip, else the test suite
    # will hang
    for test in hanging_tests:
        for tag in test['tags']:
            if (tag.compiler == compiler or tag.compiler == Compiler.all) \
                    and (tag.build == build_type or tag.build == Build_type.all):
                tests_to_skip.append(test)

    if incl_failing_tests:
        return tests_to_skip

    for test in failing_tests:
        for tag in test['tags']:
            if (tag.compiler == compiler or tag.compiler == Compiler.all) \
                    and (tag.build == build_type or tag.build == Build_type.all):
                tests_to_skip.append(test)

    return tests_to_skip


def set_tests_to_repeat(executable: str, n_repeats: int) -> dict:
    """
    Generate a dict of tests to repeat, for a given build type.

    :param str executable: executable
    :param int n_repeats: Number of times a failed job should repeat.
    :return dict tests_to_repeat: Tests to repeat, of form {name: n_repeats}.
    """
    tests_to_repeat = {}

    if n_repeats == 0:
        return tests_to_repeat

    compiler = get_compiler_type()

    if compiler is None:
        print('Compiler could not be determined from build/make.inc\n'
              'Excluding all failing tests from the test suite')
        compiler = Compiler.all

    build_type = build_type_str_to_enum[executable]

    for test in repeat_tests:
        for tag in test['tags']:
            if (tag.compiler == compiler or tag.compiler == Compiler.all) \
                    and (tag.build == build_type or tag.build == Build_type.all):
                tests_to_repeat[test['name']] = n_repeats

    return tests_to_repeat


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
