"""
Create a new test case for the test suite.

To run, for example:
  cd <EXCITINGROOT>/test
  python3 newtestcase.py -for groundstate -r path/2/test_reference -t test_farm/method/new_test

Exclude -t and the new test will be placed in the current directory:
  python3 newtestcase.py -r test_reference
"""
import sys
import os
import argparse as ap
import warnings
from typing import List

from runtest import get_test_directories
from src.exciting_settings.constants import settings
from src.io.file_system import copy_calculation_inputs
from src.reference.reference import run_single_reference
from src.reference.generate_tolerance import generate_tolerance_file
from src.runner.profile import build_type_enum_to_str, get_calculation_types, ExcitingCalculation
from src.runner.set_inputs import input_files_for_tests

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return '%s: %s \n' % (category.__name__, message)


warnings.formatwarning = warning_on_one_line


def option_parser() -> dict:
    p = ap.ArgumentParser(
        description="Usage: python3 newtestcase.py -for <method> -r <reference dir> -t <target dir>")

    p.add_argument('-for',
                   help="Method of test case, informing the choice of tolerance file."
                        "If not specified, the method attempt to be inferred from the target directory.",
                   type=str,
                   default=None)

    p.add_argument('-r',
                   metavar='--reference',
                   help="Path to an existing calculation that will be the reference.",
                   type=str,
                   default=None)

    p.add_argument('-t',
                   metavar='--target',
                   help="Target location to copy the reference to",
                   type=str,
                   default=None)

    args = p.parse_args()

    return args.__dict__


def determine_test_method(input_method: str, input_target_dir: str) -> List[ExcitingCalculation]:
    """
    Determine the calculation method.

    If the method is specified, check it's valid.
    If it is not specified, attempt to infer it from the target directory.

    :param str input_method: Input method
    :param str input_target_dir: Input target directory
    :return List[ExcitingCalculation] exciting calculation/method
    """
    if input_method is not None:
        return get_calculation_types([input_method])

    if (input_method is None) and (input_target_dir is None):
        sys.exit("Must specify either the target directory or the method")

    sub_directories = input_target_dir.split('/')
    try:
        i = sub_directories.index('test_farm')
    except ValueError:
        sys.exit("Cannot infer the method from the target directory.\n"
                 "Please specify the method at input.")

    method_directory = sub_directories[i + 1]
    return get_calculation_types([method_directory])


def set_reference_location(input_ref_location: str) -> str:
    """
    Set location of reference data that will form a new test case

    :param str input_ref_location: Input reference directory
    :return str description: Test reference directory
    """
    if input_ref_location is None:
        sys.exit("An existing calculation reference is required. "
                 "Please enter the path to the corresponding directory. \n"
                 "Make sure that a valid %s and all necessary species files are contained.\n"
                 % settings.input_file + "\n")

    reference_location = input_ref_location

    while not os.path.exists(reference_location):
        reference_location = input(
            "The reference path you entered does not exist. "
            "Please reenter or exit (type 'exit').\n")

        if reference_location.lower() == "exit":
            sys.exit()

    return reference_location


def set_target_directory(input_target_dir: str) -> str:
    """
    Set target location to copy test reference folder to

    :param str input_target_dir: Input target directory for test case
    :return str target_dir: Final target directory for test case
    """
    default_name = 'ref_data'

    if input_target_dir is None:
        print("No target has been specified to copy the reference to.")
        target_dir = os.path.join(os.getcwd(), default_name)
        print("The current directory will be used:", target_dir)
    else:
        return input_target_dir

    return target_dir


def check_clear_test_name(target_test: str):
    """
    Check if an existing test already exists with the target test name.

    target_test is some relative path, for example:
      new_test
      ./new_test
      test_farm/groundstate/new_test

    If target_test is already present in current_tests, the user will be
    given the choice to a) enter a new name, b) replace the existing test or
    c) exit.

    Note, because of test farm's subdirectory structure, it is possible to
    have two test directories with the same base name, as long as they are
    in different sub-folders of test_farm.
    """

    # Tests in test farm
    current_tests = get_test_directories(settings.test_farm, basename=False)

    # Test already exists in current operating directory
    if os.path.isdir(target_test):
        current_tests.append(target_test)

    while target_test in current_tests:
        old_target_test = target_test
        target_test = input("A test case with the same name already exists: " + old_target_test + "\n" +
                            "Please choose a different name, replace the existing test case "
                            "(enter 'replace')\n or quit (enter 'exit').\n")

        if target_test == 'replace':
            os.system('rm -r %s' % old_target_test)
            return old_target_test

        elif target_test.lower() == 'exit':
            sys.exit()

    return target_test


def main(args: dict):
    """
    Create a new test case from an input file and run the reference.

    :param args: Parsed command line arguments
    """

    calculation_types = determine_test_method(args['for'], args['t'])
    reference_location = set_reference_location(args['r'])
    target_location = set_target_directory(args['t'])
    species_files = next(os.walk(settings.species))[2]

    check_clear_test_name(target_location)
    os.makedirs(os.path.join(target_location, settings.ref_dir))

    print('Create new test case %s.' % target_location)
    print("Take reference input from %s." % reference_location)

    test_directory = os.path.join(target_location, settings.ref_dir)
    for calculation in calculation_types:
        print('Generating tolerance file for calculation type:', calculation)
        generate_tolerance_file(calculation, test_directory)

    input_files: dict = input_files_for_tests([reference_location])
    input_files_for_ref: List[str] = input_files[reference_location]
    copy_calculation_inputs(reference_location, test_directory, input_files_for_ref)
    binary = os.path.join(settings.exe_dir, build_type_enum_to_str[settings.exe_ref])
    run_single_reference(test_directory, binary, settings)


if __name__ == "__main__":
    args = option_parser()
    main(args)
