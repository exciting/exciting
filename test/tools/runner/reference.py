import os
from typing import List

from .execute import execute
from ..infrastructure import flatten_directory, remove_ignored_files

def run_single_reference(test_farm:List[str], main_out:str, test_dir:str, ref_dir:str, executable:str, 
                         ignored_output:List[str], max_time:int):
    """
    Reference run for single test case.
    :param test_farm:          location of the test farm
    :param main_out:           main output file of the exciting calculation
    :param test_dir:           test case for that the reference will be calculated
    :param ref_dir:            name of the ref directory of the test case
    :param execuatable:        executable for the exciting run
    :param ignored_output:     files that will not be saved as reference
    :param max_time:           maximum allowed time for the exciting run 
    """
    ref_path = os.path.join(test_farm,test_dir,ref_dir)

    runSuc, errMess, timing = execute(ref_path, executable, main_out, max_time)
    if not runSuc:
        print('%s: Run failed'%test_dir)
        for err in errMess:
            print(err)
        
        print('Time (s): %.1f' % timing)
        return

    print('%s: Run succeeded'%test_dir)

    # Flatten the reference directory
    flatten_directory(ref_path)

    # Append a .ref to all reference files.
    remove_ignored_files(ref_path, ignored_output)

    print('Time (s): %.1f' % timing)
    

def run_references(test_farm:str, main_out:str, testList:List[str], ref_dir:str, executable:str, 
                   ignored_output:List[str], max_time:int):
    """
    Reference run for all tests (see run_single_reference).
    :param test_farm:         location of the test farm
    :param main_out:          main output file of the exciting calculation
    :param testList:          list of test cases for that the reference will be calculated
    :param test_dir:          test case for that the reference will be calculated
    :param ref_dir:           name of the ref directory of the test case
    :param execuatable:       executable for the exciting run
    :param ignored_output:    files that will not be saved as reference
    """
    for test_dir in testList:
        run_single_reference(test_farm, 
                             main_out, 
                             test_dir, 
                             ref_dir, 
                             executable, 
                             ignored_output, 
                             max_time)



