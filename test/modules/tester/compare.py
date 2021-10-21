"""
Functions that check the type of an object
"""
from os import path
import numpy as np
from typing import Union, Tuple, List, Optional
from collections.abc import Hashable
import sys
from copy import deepcopy

from .failure import *
from .path import *
from ..termcolor_wrapper import print_color

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


def all_close_to_zero(a: np.ndarray, a_tol):
    """
    Check if an array is close to zero

    :param np.ndarray a: Numpy array
    :param a_tol: Absolute tolerance
    """
    return np.allclose(a, np.zeros(shape=a.shape), atol=a_tol)


# Comparison logic for each data type
# Each defined such that difference <= tolerance gives true
comparison_function = {int: lambda diff, tol: diff <= tol,
                       float: lambda diff, tol: diff <= tol,
                       str: lambda diff, unused_tol: diff == '',
                       list: all_close_to_zero,
                       np.ndarray: all_close_to_zero
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


class ErrorContainer:
    def __init__(self,
                 key: str,
                 test_value: Union[int, float, str, list, np.ndarray],
                 ref_value: Union[int, float, str, list, np.ndarray],
                 diff: Union[int, float, str, list, np.ndarray],
                 tol: Union[int, float, str]
                 ):
        """
        Initialise an instance of ErrorContainer class.

        For arrays, this routine only retains the errors in the array and the corresponding indices.
        np.where returns:
           indices = ((i1, i2, ...), (j1, j2, ...), ..., (n1, n2, ...))

        where (i1, j1, ..., n1) corresponds to the indices for error 1 in the test array,
        len(indices) = rank of the array, and len(indices[0]) = number of errors.

        :param str key : Concatenated key of form %key%next_key%deepest_key.
        :param Union[int, float, str, list, np.ndarray] test_value: Erroneous test data.
        :param Union[int, float, str, list, np.ndarray] ref_value: Reference data.
        :param Union[int, float, str, list, np.ndarray] diff: Difference in test and ref data.
        :param Union[int, float, str] tol: Tolerance for a difference to be considered an error.
        """
        self.key = key
        self.tol = tol

        if type(test_value) in [int, float, str]:
            self.init_with_scalars(test_value, ref_value, diff)

        elif isinstance(test_value, (list, np.ndarray)):
            # NOTE(Alex) I'm enforcing np data here as values *still* make it in as a list
            self.init_with_arrays(np.asarray(test_value), np.asarray(ref_value), np.asarray(diff))

        else:
            sys.exit("Unsupported data type for test_value in class ErrorContainer")

    def init_with_scalars(self,
                          test_value: Union[int, float, str],
                          ref_value: Union[int, float, str],
                          diff: Union[int, float, str]):
        self.test_value = test_value
        self.ref_value = ref_value
        self.diff = diff
        self.indices = None
        self.n_errors = 1

    def init_with_arrays(self,
                         test_value: np.ndarray,
                         ref_value: np.ndarray,
                         diff: np.ndarray):
        self.indices = np.where(diff > self.tol)
        self.test_value = test_value[self.indices]
        self.ref_value = ref_value[self.indices]
        self.diff = diff[self.indices]
        # Number of erroneous elements in an array
        self.n_errors = len(self.indices[0])


class ErrorFinder:
    """
    Given test data, reference data and tolerances, find any errors in the test data w.r.t. the reference data.

    The class also sorts errors by data type to make reporting (printing) informative.
    """
    def __init__(self, test_data: dict, ref_data: dict, tolerance: dict, file_name: Optional[str]=''):
        # Symbol used to indicate key concatenation
        self.cc_symbol = '%'
        assert isinstance(test_data, dict), 'test_data must be type dict'
        assert isinstance(ref_data, dict), 'test_data must be type dict'
        assert isinstance(tolerance, dict), 'test_data must be type dict'
        # List of errors, list of unused tol keys, list of keys in test_data but not in tol dict
        self.errors, self.unused_tol_keys, self.keys_not_in_tol = self.compare(test_data, ref_data, tolerance)
        # Total number of errors
        self.n_errors = len(self.errors)
        # Total number of tolerances that were not referenced
        self.n_unused_tol_keys = len(self.unused_tol_keys)
        # Total number of assertions performed
        self.n_assertions = len(tolerance) - len(self.unused_tol_keys)
        # Errors stored by data type
        self.errors_by_type = self.log_errors_by_data_type()
        # Number of errors by data type
        self.n_errors_by_type = self.count_errors_by_data_type()
        # File name
        self.file_name = file_name

    def get_dict(self) -> dict:
        """
        Convert List[ErrorContainer] into a dictionary.

        Return a dictionary of the form:
        {key1: {'test_value': test_value, 'ref_value': ref_value, 'diff': diff, 'tol': tol}, ...}

        :return dict errors_as_dict: Errors
        """
        errors_as_dict = {}
        for error in self.errors:
            error_dict = deepcopy(error.__dict__)
            key = error_dict.pop('key')
            errors_as_dict[key] = error_dict
        return errors_as_dict

    def get_error_keys(self, last=False) -> List[str]:
        """
        Get the keys associated with the errors.
        If no errors are present, None is returned.

        :param bool last: If False, return the whole key with self.cc_symbol:  '%some%nest%keys'
                          If True, return the deepest key: 'keys'
        """
        if len(self.errors) == 0:
            return None

        if last:
            get_key = lambda key: key.split(self.cc_symbol)[-1]
        else:
            get_key = lambda key: key

        return [get_key(error.key) for error in self.errors]

    def get_entry(self, key: str) -> ErrorContainer:
        """
        Given a key, get the error: ErrorContainer entry associated with it.

        :param str key: Key associated with the error.
        """
        if len(self.errors) == 0:
            return None

        error_keys = self.get_error_keys()
        try:
            index = error_keys.index(key)
            return self.errors[index]
        except KeyError:
            raise KeyError('Key not found in errors: ' + key)

    def get_error_value(self, key: str) -> Union[int, float, str, np.ndarray]:
        """
        Given a key, get the error value entry associated with it.

        :param str key: Key associated with the error.
        """
        if len(self.errors) == 0:
            return None
        error_obj = self.get_entry(key)
        return error_obj.diff

    def log_errors_by_data_type(self) -> dict:
        """
        Logs the keys of regression errors by data type.

        :return dict error_by_type: Dict containing error keys by data type.
        """
        errors_by_type = {int: [], float: [], str: [], list: [], np.ndarray: []}

        if self.n_errors == 0:
            return errors_by_type

        for error in self.errors:
            data_type = type(error.test_value)
            errors_by_type[data_type].append(error)

        return errors_by_type

    def count_errors_by_data_type(self) -> dict:
        """
        Count errors by their data type.

        :return dict n_errors_by_type: Dict containing number of errors by data type.
        """
        assert len(self.errors_by_type[list]) == 0, "All lists should be converted to np array"

        n_errors_by_type = {'integer': len(self.errors_by_type[int]),
                            'float': len(self.errors_by_type[float]),
                            'string': len(self.errors_by_type[str]),
                            'array':  len(self.errors_by_type[np.ndarray]),
                            }

        total = sum([n for n in n_errors_by_type.values()])
        n_errors_by_type['total'] = total

        return n_errors_by_type

    def compare(self, test_data: dict, ref_data: dict, tolerance: dict) -> Tuple[List[ErrorContainer], List[str], List[str]]:
        """
        Wrapper function for recursive comparison of a dictionary of test data with a dictionary of reference data.

        Valid test_data and ref_data dictionary values are:
        int, float, str, List[int], List[float], List[dict], np.ndarray and dict.

        For a full list of hashable data types, and containers, see:
        https://stackoverflow.com/questions/14535730/what-does-hashable-mean-in-python

        :param dict test_data: Test data.
        :param dict ref_data: Reference data, with same keys and nesting as test_data.
        :param dict tolerance: Tolerances for keys with hashable values.

        :return List[ErrorContainer] errors: Key:values of test values that exceeded reference + tolerance values.
                             If no errors are found, the routine returns {}.
                             For comparison of lists and arrays, if any errors are found, all the diffs are returned.
        :return List[str] unused_tol_keys: Keys of any tolerances that are not evaluated.
        :return List[str] keys_not_in_tolerance: Keys of any test data not evaluated.
        """
        errors: List[ErrorContainer] = []
        used_tol_keys: set = set()
        keys_not_in_tolerance: set = set()

        errors, used_tolerances, keys_not_in_tolerance = \
            self._recursive_compare_reference_with_target(test_data,
                                                          ref_data,
                                                          tolerance,
                                                          errors,
                                                          used_tol_keys,
                                                          keys_not_in_tolerance,
                                                          full_key='')
        all_tol_keys = set(tolerance.keys())

        return errors, list(all_tol_keys - used_tol_keys), list(keys_not_in_tolerance)

    def _recursive_compare_reference_with_target(self,
                                                 test_data,
                                                 ref_data,
                                                 tolerance: dict,
                                                 errors: list,
                                                 used_tol_keys: set,
                                                 keys_not_in_tolerance: set,
                                                 full_key='') -> Tuple[list, set, set]:
        """
        Recursively compare a dictionary of test data with a dictionary of reference data.
        This routine should is not intended to be called outside of this module.

        Description
        --------------
        Passes test_data and ref_data recursively, but the whole tolerance dict is always
        passed as it should only contain keys for hashable values.

        Valid test_data and ref_data dictionary values are:
        int, float, str, List[int], List[float], List[dict] and dict.

        For a full list of hashable data types, and containers, see:
        https://stackoverflow.com/questions/14535730/what-does-hashable-mean-in-python

        Notes
        --------------
        Convoluted Data Structures
        Poorly-considered data structures should not be facilitated by this function,
        implying exciting's parsers should be sensible. For example one could add:

          {..., 'd':[1, 2, {'e':3, 'f': 4}]}

        to this function's logic, but there shouldn't be a reason for a parser to store data in this manner.

        TODO(A/B/H) Issue 97. Unevaluated Test Data Keys
        Deliberately Ignored Test Data Keys vs Erroneously Ignored
        Some parsed data are deliberately not evaluated. For examples, 'Species', 'parameters loaded from' and
        'Wall time (seconds)' (amongst others) in INFO.OUT.
        raise KeyError only wants to raise an error if a key has accidentally been missed out of the tolerance file.

        A nicer approach is to preprocess the keys and remove ones that should not be evaluated from the
        test (and reference) data prior to passing the to this routine. Then keyError can be raised and
        remove `keys_not_in_tolerance`

        Arguments
        --------------
        :param test_data: Test data.
        :param ref_data: Reference data.
        :param dict tolerance: Tolerances for keys with hashable values.
        :param dict errors: Errors: Keys of any tolerances that are not evaluated.
        :param set used_tol_keys: Keys of any test data not evaluated.
        :param str full_key: Concatenated str of keys, for a given error.

        :return list errors: List of ErrorContainer objects.
                             If no errors are found, the routine returns [].
        :return set used_tol_keys: Keys of any tolerances that are not evaluated.
        :return set keys_not_in_tolerance: Keys of any test data not evaluated.
        """
        assert type(test_data) == type(ref_data), "test_data and ref_data are different types"
        assert len(test_data) == len(ref_data), "Length of test_data differs from length of ref_data"

        for key, test_value in test_data.items():
            full_key += self.cc_symbol + str(key)

            if isinstance(test_value, dict):
                self._recursive_compare_reference_with_target(test_value, ref_data[key], tolerance, errors,
                                                              used_tol_keys, keys_not_in_tolerance, full_key)

            if isinstance(test_value, list) and (not all_hashable_or_dict(test_value)):
                # NOTE (Alex) I added this because for the GW parser, data parsed as np.ndarray was
                # ending up as a list. This circumvents the problem but doesn't fix the cause (and I don't have
                # the time or insight to debug)
                try:
                    test_value = np.asarray(test_value)
                except ValueError:
                    raise ValueError('All elements of a parsed list should be hashable or dict')

            if isinstance(test_value, (int, float, str, np.ndarray)) or hashable_list(test_value):
                difference_function = diff_condition[type(test_value)]
                diff = difference_function(test_value, ref_data[key])
                all_close = comparison_function[type(test_value)]

                try:
                    tol = tolerance[key]
                    # Log all evaluated tolerance keys
                    used_tol_keys.add(key)

                    if not all_close(diff, tol):
                        errors.append(ErrorContainer(full_key.lstrip(self.cc_symbol), test_value, ref_data[key], diff, tol))

                except KeyError:
                    # Key not present in tolerance file. This implies either:
                    #   a) Data deliberately not checked.
                    #   b) New data has been added to an exciting output and parsed but its key has
                    #      not been added to the tolerance file.
                    keys_not_in_tolerance.add(key)

            full_key = self.remove_last_key_string(full_key)

        return errors, used_tol_keys, keys_not_in_tolerance

    def remove_last_key_string(self, full_key: str):
        """
        For a key %a%b%c%, roll back to %a%b.

        Used when full_key has been summed to, but the last concatenation was
        of a key that is on the same level of nesting.

        Example:
          full_key = "%a%b%c%"
          full_key.split(self.cc_symbol)[1:-1] = [a, b, c]
           "".join(self.cc_symbol + s for s in [a, b, c]][:-1]) = "%a%b"
          return  "%a%b"
        """
        return "".join(self.cc_symbol + s for s in full_key.split(self.cc_symbol)[1:-1])

    def print_scalar_errors(self):
        """
        Print scalar (single value) errors, where scalars are defined as int, float and string.
        """
        n_scalar_errors = sum([self.n_errors_by_type[key] for key in ['integer', 'float', 'string']])
        if n_scalar_errors == 0:
            return

        indent = ' ' * 6
        failure_msg = {int: 'INTEGER FAILURES', float: 'FLOAT FAILURES', str: 'STRING FAILURES'}
        header_format = {int: '{:<40} {:>13} {:>13} {:>13} {:>13}', float: '{:<40} {:>13} {:>13} {:>13} {:>13}', str: '{:<40} {:>34} {:>34}'}
        data_format = {int: '{:<40} {:>13} {:>13} {:>13} {:>13}', float: '{:<40} {:13.8f} {:13.8f} {:13.2e} {:13.2e}', str: '{:<40} {:>34} {:>34}'}

        for data_type in [int, float,str]:
            errors = self.errors_by_type[data_type]
            if len(errors) == 0:
                continue

            print(indent + failure_msg[data_type])
            print(indent + header_format[data_type].format('Key', 'Test Data', 'Ref Data', 'Diff', 'Tolerance'))
            print(indent + '-' * 110)

            for error in errors:
                format_specifier = data_format[data_type]
                print(indent + format_specifier.format(error.key, error.test_value, error.ref_value, error.diff, error.tol))
            print()

    def print_array_errors(self):
        """
        Print array errors, where array is defined as `np.ndarray`, only.
        """
        if self.n_errors_by_type['array'] == 0:
            return

        indent = ' ' * 6
        header_format = '{:<40} {:>13} {:>13} {:>13} {:>13}'
        data_format = '{:<40} {:13.8f} {:13.8f} {:13.2e} {:13.2e}'

        label = {1: '(element)',
                 2: '(row, column)',
                 3: '(i, j, k)',
                 4: '(i, j, k, l)',
                 5: '(i, j, k, l, m)',
                 6: '(i, j, k, l, m, n)'
                 }

        print('      ARRAY FAILURES')

        for array in self.errors_by_type[np.ndarray]:
            dimensions = len(array.test_value.shape)
            try:
                indices_label = label[dimensions]
            except KeyError:
                raise KeyError('Array dimension exceeds 6')

            print(indent + 'Key:', array.key)
            print(indent + header_format.format(indices_label, 'Test Data', 'Ref Data', 'Diff', 'Tolerance'))
            print(indent + '-' * 110)

            # Print all erroneous elements for a given test array
            for i in range(0, array.n_errors):
                indices_str = "(" + ",".join(str(element[i]) for element in array.indices) + ")"
                print(indent + data_format.format(indices_str, array.test_value[i], array.ref_value[i], array.diff[i], array.tol))

            print('')

    def print_result(self):
        """
        Print results of a regression test on self.file_name

        For example, a string of the format:
        ```
            info.xml SUCCESS

        ```
        for successful comparison, or
        ```
           Test Case: test_functions/dummy_app_tests/LDA_VWN-He_failing_scalars
           Time (s): 0.0
           Test execution completed
           Output Files: SUCCESS 4/5, NOT EVALUATED 0/5, FAIL 1/5.
               INFO.OUT FAIL
                 Failures TOTAL: 3, INTEGER: 1, FLOAT: 1, ARRAY: 0, STRING 1
        ```
        when failure/s occit
        """
        if self.n_errors == 0:
            print_color('    ' + self.file_name + ' SUCCESS', 'green')
            return

        description = (self.n_errors_by_type['total'],
                       self.n_errors_by_type['integer'],
                       self.n_errors_by_type['float'],
                       self.n_errors_by_type['array'],
                       self.n_errors_by_type['string'])

        print_color('    ' + self.file_name + ' FAIL', 'red')
        print_color('      Failures TOTAL: %i, INTEGER: %i, FLOAT: %i, ARRAY: %i, STRING %i' % description, 'red')

        self.print_scalar_errors()
        self.print_array_errors()





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
