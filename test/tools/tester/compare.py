"""
Functions that check the type of an object
"""
from os import path
import numpy as np
from typing import Union, List
from collections import Hashable
import sys

from .failure import *
from .path import *


#########################
# Code to retain
#########################

def strings_equal(x: str, y: str, error_mgs='strings differ', ignore_lr_whitespace=True) -> str:
    """
    Return error message if two strings differ

    :param str x: Target string
    :param str y: Reference string
    :param str error_mgs: Error message if x and y differ
    :param bool ignore_lr_whitespace: Ignore left and right whitespace when comparing strings x and y.
    """
    if ignore_lr_whitespace:
        diff = x.strip() != y.strip()
    else:
        diff = x != y

    if diff:
        return error_mgs
    return ''

# Definition of the difference in two values, for a given data type
diff_condition = {int: lambda x, y: abs(x - y),
                  float: lambda x, y: abs(x - y),
                  str: strings_equal,
                  list: lambda x, y: np.abs(np.array(x) - np.array(y)),
                  np.ndarray: lambda x, y: np.abs(x - y)
                  }


# Comparison logic for each data type
# Each defined such that difference <= tolerance gives true
comparison_function = {int: lambda diff, tol: diff <= tol,
                       float: lambda diff, tol: diff <= tol,
                       str: lambda diff, unused_tol: diff == '',
                       list: np.allclose,
                       np.ndarray: np.allclose
                       }


def all_hashable(data) -> bool:
    """
    Checks if all elements of data are hashable
    :param data: Container
    """
    return all(isinstance(x, Hashable) for x in data)


def hashable_list(my_list: list) -> bool:
    """
    Checks if all elements of data are hashable
    :param list my_list: List to check
    """
    return isinstance(my_list, list) and all_hashable(my_list)


def all_hashable_or_dict(my_list: list) -> bool:
    """
    Checks if all elements in a list are hashable or a dictionary.
    :param list my_list: List to check
    """
    return all(isinstance(x, (Hashable, dict)) for x in my_list)


def compare_reference_with_target(test_data: dict, ref_data: dict, tolerance: dict) -> Union[dict, List[str]]:
    """
    Wrapper function for recursive comparison of a dictionary of test data with a dictionary of reference data.

    Valid test_data and ref_data dictionary values are:
    int, float, str, List[int], List[float], List[dict] and dict.

    For a full list of hashable data types, and containers, see:
    https://stackoverflow.com/questions/14535730/what-does-hashable-mean-in-python

    :param dict test_data: Test data.
    :param dict ref_data: Reference data, with same keys and nesting as test_data.
    :param dict tolerance: Tolerances for keys with hashable values.

    :return dict errors: Key:values of test values that exceeded reference + tolerance values.
                         If no errors are found, the routine returns {}.
                         For comparison of lists and arrays, if any errors are found, all the diffs are returned.
    :return List[str] unused_tol_keys: Keys of any tolerances that are not evaluated.
    """

    # Implies an error in the caller
    assert type(test_data) == type(ref_data), "test_data and ref_data are different types"

    # Implies a change in the output data
    assert len(test_data) == len(ref_data), "Length of test_data differs from length of ref_data"

    errors = {}
    used_tol_keys = set()
    errors, used_tolerances = _compare_reference_with_target_implementation(test_data,
                                                                            ref_data,
                                                                            tolerance,
                                                                            errors,
                                                                            used_tol_keys)
    all_tol_keys = set(tolerance.keys())

    return errors, list(all_tol_keys - used_tol_keys)


# Private function. Should only be called via compare_reference_with_target (unless one knows that they've doing)
def _compare_reference_with_target_implementation(test_data,
                                                  ref_data,
                                                  tolerance: dict,
                                                  errors: dict,
                                                  used_tol_keys: set) -> Union[dict, set]:
    """
    Recursively compare a dictionary of test data with a dictionary of reference data.
    This routine should is not intended to be called outside of this module.

    Passes test_data and ref_data recursively, but the whole tolerance dict is always
    passed as it should only contain keys for hashable values.

    Valid test_data and ref_data dictionary values are:
    int, float, str, List[int], List[float], List[dict] and dict.

    For a full list of hashable data types, and containers, see:
    https://stackoverflow.com/questions/14535730/what-does-hashable-mean-in-python

    Note: Poorly-considered data structures should not be facilitated by this function,
    implying exciting's parsers should be sensible. For example one could add:

      {..., 'd':[1, 2, {'e':3, 'f': 4}]}

    to this function's logic, but there shouldn't be a reason for a parser to store data
    in this manner.

    :param test_data: Test data.
    :param ref_data: Reference data.
    :param dict tolerance: Tolerances for keys with hashable values.
    :param dict errors: Errors
    :param set used_tol_keys:

    :return dict errors: Key:values of test values that exceeded reference + tolerance values.
                         If no errors are found, the routine returns {}.
                         For comparison of lists and arrays, if any errors are found, all the diffs are returned.
    :return set used_tol_keys: Keys of any tolerances that are not evaluated.
    """

    for key, value in test_data.items():
        if isinstance(value, dict):
            _compare_reference_with_target_implementation(value, ref_data[key], tolerance, errors, used_tol_keys)

        if isinstance(value, list) and (not all_hashable_or_dict(value)):
            raise ValueError('All elements of a parsed list should be hashable or dict')

        if isinstance(value, (int, float, str, np.ndarray)) or hashable_list(value):
            difference_function = diff_condition[type(value)]
            diff = difference_function(value, ref_data[key])
            all_close = comparison_function[type(value)]

            if not all_close(diff, tolerance[key]):
                errors[key] = diff

            # Log all evaluated tolerance keys
            used_tol_keys.add(key)

    return errors, used_tol_keys


#################################################################
# Code to scrap once JSON tols are fully implemented
#################################################################


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
