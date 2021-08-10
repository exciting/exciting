"""
Functions that check the type of an object
"""
from os import path
import numpy as np
from typing import Union
from typing import List

from .failure import *
from .path import *



error_condition = {int: lambda x, y: abs(int(x) - int(y)),
                   float: lambda x, y: abs(float(x) - float(y)),
                   list: lambda x, y: np.sqrt(np.mean((np.array(x) - np.array(y)) ** 2)),
                   np.ndarray: lambda x, y: np.sqrt(np.mean((np.array(x) - np.array(y)) ** 2))
                   }


def set_tolerance(test_failures: List[Failure], ref_data):
    """
    Chooses appropriate tolerance.
    :param list test_failures:  List of test failures.
    :param ref_data:            Reference output data from test
    """
    if test_failures.path.lastElement() in test_failures.tolValuePair['value']:
        return(float(test_failures.tolValuePair['tol']))
    else:
        default_tolerance = {int: test_failures.tolFloat,
                             float: test_failures.tolFloat,
                             list: test_failures.tolMSE,
                             np.ndarray: test_failures.tolMSE,
                             str: None
                             }
        return(default_tolerance[type(ref_data)])


def check_fail(test_failures: List[Failure], 
               test_data, 
               ref_data, 
               error_condition: dict, 
               test_dir:path, 
               test_name: str, 
               tol = None, 
               path = None):
    """
    Check if a path (key) agrees between test_data and ref_data, given a tolerance.

    :param list test_failures:      List of test failures
    :param test_data:               Output data from test
    :param ref_data:                Reference output data from test
    :param dict error_condition:    Given type as a key, return the appropriate error function
    :param path test_dir:           Path to test case directory
    :param str test_name:           Name of output file
    """
    if tol == None:
        tol = set_tolerance(test_failures, ref_data)
    
    if path == None:
        path = test_failures.path

    if type(test_data) == str:
        check = test_data != ref_data
        error = "Strings differ"
    else:
        check = error_condition[type(test_data)](test_data, ref_data) > tol
        error = error_condition[type(test_data)](test_data, ref_data)

    if check:
        test_failures.append(Failure(error=error,
                                     tolerance=tol,
                                     path=path,
                                     test_data=test_data,
                                     ref_data=ref_data,
                                     test_dir = test_dir,
                                     test_name = test_name
                                     )
                            )
    return test_failures
   

def compare_data(test_data: Union[dict, list], 
                 ref_data: Union[dict, list], 
                 test_dir: str, 
                 test_failures: List[Failure]):
    """
    :param dict/list test_data: Output data from test
    :param dict/list ref_data:  Reference output data from test
    :param str test_dir:        Path to test directory
    :param list test_failures:           List of test failures
    """

    if isinstance(test_data, list) and isinstance(ref_data, list):
        if len(test_data) != len(ref_data):
            test_failures.append(Failure(failure_code=Failure_code.FORMAT, path=test_failures.path))
            return

        all_floats_in = lambda data: all(isinstance(x, float) for x in data)
        if all_floats_in(test_data) and all_floats_in(ref_data):
            check_fail(test_failures, test_data, ref_data, error_condition, test_dir, test_failures.name)
        else:
            for i in range(0, len(test_data)):
                test_failures.path = test_failures.path.update(str(i))
                if isinstance(test_data[i], dict) and isinstance(ref_data[i], dict):
                    compare_data(test_data[i], ref_data[i], test_dir, test_failures)

                elif isinstance(test_data[i], list) and isinstance(ref_data[i], list):
                    compare_data(test_data[i], ref_data[i], test_dir, test_failures)

                elif type(test_data[i])==type(ref_data[i]):
                    check_fail(test_failures, test_data[i], ref_data[i], error_condition, test_dir, test_failures.name)

                else:
                    test_failures.append(Failure(failure_code=Failure_code.FORMAT, path=test_failures.path))
                    return
                test_failures.path = test_failures.path.removeLastElement()

    elif isinstance(test_data, dict) and isinstance(ref_data, dict):
        if test_data.keys() != ref_data.keys():
            test_failures.append(Failure(failure_code=Failure_code.FORMAT, path=test_failures.path))
        for key in test_data.keys():
            if key not in ref_data.keys():
                test_failures.ignore.append(key)
        for key in ref_data.keys():
            if key not in test_data.keys():
                 test_failures.ignore.append(key)
        for key in ref_data.keys():
            if key in test_failures.ignore:
                break
            assert isinstance(key, (str, int, float)), "Result key not suitable for convertion to string"
            test_failures.path = test_failures.path.update(str(key))
            if isinstance(test_data[key], dict) and isinstance(ref_data[key], dict):
                compare_data(test_data[key], ref_data[key], test_dir, test_failures)

            elif isinstance(test_data[key], list) and isinstance(ref_data[key], list):
                compare_data(test_data[key], ref_data[key], test_dir, test_failures)

            elif type(ref_data[key])==type(test_data[key]):
                check_fail(test_failures, test_data[key], ref_data[key], error_condition, test_dir, test_failures.name)

            else:
                test_failures.append(Failure(failure_code=Failure_code.FORMAT, path=test_failures.path))
                return
            test_failures.path = test_failures.path.removeLastElement()


    else:
        test_failures.append(Failure(failure_code=Failure_code.FORMAT, path=test_failures.path))


def compare_INFO(test_data: Union[dict, list], 
                 ref_data: Union[dict, list],
                 test_dir: str,
                 test_failures: List[Failure],
                 tolIter: int,
                 maxIter: int):
    """
    :param dict/list test_data: Output data from test
    :param dict/list ref_data:  Reference output data from test
    :param str test_dir:        Path to test directory
    :param list test_failures:  List of test failures
    :param int tolIter:         Tolerance for number of iterations
    :param int tolIter:         Maximum number of iterations
    """

    test_max_scf = find_max_scf(test_data['scl'])
    ref_max_scf = find_max_scf(ref_data['scl'])

    # checks whether the number of iterations is in the tolerance
    if maxIter != 0:
        if test_max_scf > int(maxIter):
            test_failures.append(
                Failure(error=max(test_max_scf) - int(maxIter), tolerance=0, path="INFO.OUT/Iterations"))
    else:
        check_fail(test_failures, test_max_scf, ref_max_scf, error_condition, test_dir, test_failures.name, tol=int(tolIter), path="INFO.OUT/Iterations")


    # TODO(A/B/H)issue 67. Add integer failure

    # compares the last iteration
    compare_data(test_data['scl'][str(test_max_scf)], ref_data['scl'][str(ref_max_scf)], test_dir, test_failures)

    # compares the initialization
    compare_data(test_data['initialization'], ref_data['initialization'], test_dir, test_failures)


def compare_info(test_data: Union[dict, list], 
                 ref_data: Union[dict, list],
                 test_dir: str,
                 test_failures: List[Failure],
                 tolIter: int,
                 maxIter: int):
    """
    :param dict/list test_data: Output data from test
    :param dict/list ref_data:  Reference output data from test
    :param str test_dir:        Path to test directory
    :param list test_failures:  List of test failures
    :param int tolIter:         Tolerance for number of iterations
    :param int maxIter:         Maximum number of iterations
    """
    
    test_max_scf = find_max_scf(test_data['scl'])
    ref_max_scf = find_max_scf(ref_data['scl'])

    # checks whether the number of iterations is in the tolerance
    if maxIter != 0:
        if test_max_scf > int(maxIter):
            test_failures.append(
                Failure(error=test_max_scf - int(maxIter), tolerance=0, path="info.xml/Iterations"))
    else:
        check_fail(test_failures, test_max_scf, ref_max_scf, error_condition, test_dir, test_failures.name, tol=int(tolIter), path="info.xml/Iterations")

    # compares the energies calculated in the last iteration
    compare_data(test_data['scl'][str(test_max_scf)]['energies'], ref_data['scl'][str(ref_max_scf)]['energies'], test_dir,
                 test_failures)

    # if moments are calculated, they are compared too
    if 'moments' in ref_data['scl'][str(ref_max_scf)]:
        compare_data(test_data['scl'][str(test_max_scf)]['moments'], ref_data['scl'][str(ref_max_scf)]['moments'], test_dir,
                     test_failures)


def compare_eigval(test_data: Union[dict, list], 
                   ref_data: Union[dict, list], 
                   test_dir: str, 
                   test_failures: List[Failure],   
                   eigval: dict):
    """
    :param dict/list test_data: Output data from test
    :param dict/list ref_data:  Reference output data from test
    :param str test_dir:        Path to test directory
    :param list test_failures:  List of test failures
    :param dict eigval:         Dictionary containing pairs of eigenvalues
    """
    compare_data(test_data, ref_data, test_dir, test_failures)
    if 'eigvalPair' in eigval:
        for key in eigval['eigvalPair'].split():
            key1 = key.split(',')
            ref = float(ref_data['kpt']['1']['state'][key1[1]]['eigenvalue']) - float(
                ref_data['kpt']['1']['state'][key1[0]]['eigenvalue'])
            dat = float(test_data['kpt']['1']['state'][key1[1]]['eigenvalue']) - float(
                test_data['kpt']['1']['state'][key1[0]]['eigenvalue'])
            path1 = Path('eigval.xml/kpt/1/state/' + str(key1[1]) + '-' + str(key1[0]) + '/eigenvalue_diff')
            check_fail(test_failures, dat, ref, error_condition, test_dir, test_failures.name, tol=float(eigval['tolEigval']), path=str(path1))


def find_max_scf(scf_results: dict) -> int:
    """
    Returns the number of scf loops 
    :param dict scf_results: dictionary containing all scf loops
    :return: int:            total number of scf loops
    """
    scf_keys = [int(key) for key in scf_results.keys()]
    return max(scf_keys)
