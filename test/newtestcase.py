"""
Create a new test case for the test suite.

To run, for example:
  cd <EXCITINGROOT>/test
  python3 newtestcase.py -r path/2/test_reference -t test_farm/groundstate/new_test -i init_groundstate.xml

Exclude -t and the new test will be placed in the current directory:
  python3 newtestcase.py -r test_reference -i init_groundstate.xml
"""
import sys
import os
import argparse as ap
import warnings

from runtest import get_test_directories
from modules.constants import settings
from modules.infrastructure import copy_exciting_input, create_init
from modules.runner.reference import run_single_reference
from modules.utils import build_type_enum_to_str


def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return '%s: %s \n' % (category.__name__, message)


warnings.formatwarning = warning_on_one_line


def option_parser() -> dict:
    p = ap.ArgumentParser(
        description="Usage: python3 newtestcase.py -d <test description> -r <reference dir> -t <target dir> "
                    "-i <init_template.xml>")

    p.add_argument('-d',
                   metavar='--description',
                   help="Description of the new test case.",
                   type=str,
                   default='')

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

    p.add_argument('-i',
                   metavar='--init_file',
                   help="Init file for the new test case. "
                        "The available init files can be found at 'xml/init_templates/'. "
                        "Default is init_groundstate.xml.",
                   type=str,
                   default='init_groundstate.xml')

    args = p.parse_args()

    return args.__dict__


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

    description = args['d']
    reference_location = set_reference_location(args['r'])
    target_location = set_target_directory(args['t'])
    init_file = args['i']
    species_files = next(os.walk(settings.species))[2]

    check_clear_test_name(target_location)
    os.makedirs(os.path.join(target_location, settings.ref_dir))
    create_init(target_location, description, init_file)

    print('Create new test case %s.' % target_location)
    print("Take reference input from %s." % reference_location)

    test_directory = os.path.join(target_location, settings.ref_dir)
    copy_exciting_input(reference_location, test_directory, species_files, settings.input_file)
    binary = os.path.join(settings.exe_dir, build_type_enum_to_str[settings.exe_ref])
    run_single_reference(test_directory, binary, settings)


if __name__ == "__main__":
    args = option_parser()
    main(args)
