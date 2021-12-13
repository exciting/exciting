"""
Set input files
"""
import os
from typing import List

from ..exciting_settings.constants import settings, species_files as all_species_files


def test_case_species_files(tests: List[str], subdirectory: str) -> dict:
    """
    List the species files associated with each test case.
    Note, one could alternatively extract the information from input.xml

    :param List[str] tests: Test cases.
    :param str subdirectory: Subdirectory of test case, in which to look for input files.
    :return dict species_per_test: Species files present for each test case.
    """
    species_per_test = {}
    for test in tests:
        directory = os.path.join(test, subdirectory)
        try:
            files_in_directory = next(os.walk(directory))[2]
        except StopIteration:
            raise StopIteration(f'Directory does not exist or is empty: {directory}')

        test_case_species = list(set(files_in_directory) & set(all_species_files))
        if not test_case_species:
            raise FileNotFoundError(f"Test case reference data is missing species files: {test}")

        species_per_test[test] = test_case_species

    return species_per_test


def input_files_for_tests(tests: List[str], subdirectory='') -> dict:
    """
    Set input files for test cases

    Sets input.xml and species files.

    When restarting tests from STATE.OUT becomes possible in the suite,
    this will need extending, and probably specifying in a config file.
    """
    species_files_per_test = test_case_species_files(tests, subdirectory)
    input_xml = settings.input_file
    input_files = {test: species_files + [input_xml] for test, species_files in species_files_per_test.items()}
    return input_files
