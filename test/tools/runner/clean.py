import shutil
import os

from ..infrastructure import remove_ignored_files

def clean_tests(test_farm:str, test_list:list, run_dir:str, ref_dir:str, ignored_output:list):
    """
    Cleans the tests in test_list (see clean_single_test).
    :param test_farm:          location of the test farm
    :param test_list:          test cases that will be cleaned
    :param run_dir:            name of the run directory of the test case
    :param ref_dir:            name of the reference directory
    :param ignored_output:     Files that are ignored by the tests.
    """

    for test_dir in test_list:
        # remove run directory
        try:
            shutil.rmtree(os.path.join(test_farm, test_dir, run_dir))
        except Exception:
            pass

        # clean up ref directory
        remove_ignored_files(os.path.join(test_farm, test_dir, ref_dir), ignored_output)
    
    print('Test directories cleaned.')

