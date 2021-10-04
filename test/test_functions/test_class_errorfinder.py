"""
Test the class ErrorFinder

To run with pytest:
```
cd test
pytest -s selftests/test_class_errorfinder.py
```
"""

import unittest
import numpy as np

from ..modules.tester.compare import hashable_list, ErrorFinder


def round_values(input_data: dict, dp=14) -> dict:
    """
    Helper function to round errors to some specified decimal place.
    For example  3.9999....8 should be 4.0, but isn't due to floating point error.
    Required as unittest does not have a method assertListAlmostEqual

    Less checks and balances as this should be guaranteed by the compare function.

    :param dict input_data: Input data, with values to round
    :param int dp: Number of decimal places to round to
    :param dict input_data: Input data values are mutated - rounded to 14 dp
    """
    for key, value in input_data.items():

        if isinstance(value, dict):
            round_values(input_data)

        if isinstance(value, (int, float, np.ndarray)) or hashable_list(value):
            input_data[key] = np.round(value, dp)

    return input_data


class TestErrorFinderReferenceWithTarget(unittest.TestCase):
    """
    Given some target data in a dictionary structure, compare it to reference data
    in the same structure, and compare each difference in hashable values to a tolerance.
    """

    def test_test_and_ref_data_equal(self):
        """
        Expect each test to pass as deliberately specifying reference = tolerance
        in the argument calls
        """
        # Each hashable data type
        data = {'a': 1, 'b': 1.1, 'c': 'hello'}
        tolerance = {'a': 0, 'b': 1.e-10, 'c': ''}
        result = ErrorFinder(data, data, tolerance)
        self.assertEqual(result.errors, [], msg='Each hashable data type in compare_reference_with_target')
        self.assertEqual(result.unused_tol_keys, [], "Expected all tolerances to be checked")
        self.assertEqual(result.keys_not_in_tol, [], "Expect all keys in test data to also be in tolerances")

        # Value = string, with lr whitespace ignored
        data = {'a': '  hello  '}
        tol = {'a': ''}
        ref = {'a': 'hello'}
        result = ErrorFinder(data, ref, tol)
        self.assertEqual(result.errors, [], msg="Left and right whitespace ignored in string comparisons")
        self.assertEqual(result.unused_tol_keys, [], "Expected all tolerances to be checked")
        self.assertEqual(result.keys_not_in_tol, [], "Expect all keys in test data to also be in tolerances")

        # Value = list
        data = {'a': 1.0, 'b': [1.0, 2.0]}
        tolerance = {'a': 1.e-10, 'b': 1.e-10}
        result = ErrorFinder(data, data, tolerance)
        self.assertEqual(result.errors, [], msg="Value in target dictionary = list")
        self.assertEqual(result.unused_tol_keys, [], "Expected all tolerances to be checked")
        self.assertEqual(result.keys_not_in_tol, [], "Expect all keys in test data to also be in tolerances")

        # Value = dict
        data = {'a': 1.0, 'b': {'c': 1, 'd': 1.1, 'e': 'hello'}}
        tolerance = {'a': 1.e-10, 'c': 0, 'd': 1.e-10, 'e': ''}
        result = ErrorFinder(data, data, tolerance)
        self.assertEqual(result.errors, [], msg="Value in target dictionary = dictionary")
        self.assertEqual(result.unused_tol_keys, [], "Expected all tolerances to be checked")
        self.assertEqual(result.keys_not_in_tol, [], "Expect all keys in test data to also be in tolerances")

        # Value = Nested list and dict
        data = {'a': {'b': 1, 'c': 1.1, 'd': 'hello', 'e': {'f': 5.5}, 'g': [7, 8.9, 1.e-9]}}
        tolerance = {'b': 0, 'c': 1.e-10, 'd': '', 'f': 1.e-10, 'g': 1.e-10}
        result= ErrorFinder(data, data, tolerance)
        self.assertEqual(result.errors, [], msg="Values in target dictionary contain a mix of dict and list")
        self.assertEqual(result.unused_tol_keys, [], "Expected all tolerances to be checked")
        self.assertEqual(result.keys_not_in_tol, [], "Expect all keys in test data to also be in tolerances")

    def test_ignored_tolerance_keys(self):
        data = {'a': 1, 'b': 1.1}
        tolerance = {'a': 0, 'b': 1.e-10, 'c': np.infty}
        result = ErrorFinder(data, data, tolerance)
        self.assertEqual(result.unused_tol_keys, ['c'], msg="Keys of unevaluated tolerances")
        self.assertEqual(result.errors, [], msg="Expect errors to be empty")
        self.assertEqual(result.keys_not_in_tol, [], "Expect all key in data to also be in tolerance")

    def test_string_tolerances(self):
        data = {'a': 'Some string'}
        tolerance = {'a': 'arbitrary string tolerance'}
        result = ErrorFinder(data, data, tolerance)
        self.assertEqual(result.unused_tol_keys, [], msg="Expected key with string tolerance to be accessed")
        self.assertEqual(result.errors, [], msg="String tolerance values are not evaluated. Only the keys are used")
        self.assertEqual(result.keys_not_in_tol, [], "Expect all keys in data to also be in tolerance, hence this should be empty")

    def test_when_hashable_values_do_not_match_ref(self):
        """
        Test when the test/target values do not match the reference, to within tolerance.
        int, float and str
        """
        # Each hashable data type
        data = {'a': 1, 'b': 1.1, 'c': 'hello'}
        ref = {'a': 2, 'b': 1.2, 'c': 'goodbye'}
        tolerance = {'a': 0, 'b': 1.e-10, 'c': ''}
        result = ErrorFinder(data, ref, tolerance)

        # Extract the errors from the result object
        errors = {key.replace('%', ''): value['diff'] for key, value in result.get_dict().items()}
        errors = round_values(errors, dp=14)

        self.assertEqual(errors, {'a': 1, 'b': 0.1, 'c': 'strings differ'}, msg="Expect errors in keys a, b and c")
        self.assertEqual(result.unused_tol_keys, [], "Expected all tolerances to be checked")
        self.assertEqual(result.keys_not_in_tol, [], "Expect all key in data to also be in tolerance")

    def test_when_whole_list_does_not_match_ref(self):
        """
        Test when the test/target values do not match the reference, to within tolerance.
        list
        """
        # All list values in data differ to reference
        data = {'a': 1.0, 'b': [1.0, 2.0]}
        ref = {'a': 1.0, 'b': [1.1, 2.5]}
        tolerance = {'a': 1.e-10, 'b': 1.e-10}

        result = ErrorFinder(data, ref, tolerance)

        # Assert the key of the error
        self.assertEqual(result.get_error_keys(), ['b'], msg="Only expect error in key b")
        self.assertEqual(result.get_error_keys(last=True), ['b'], msg="No nesting of keys, so last=True returns same result")

        # Assert error value
        error = result.get_error_value('b')
        error = np.round(error, 14)
        self.assertListEqual(list(error), [0.1, 0.5], msg="All values in list differ to ref")

        # Assert other attributes
        self.assertEqual(result.unused_tol_keys, [], "Expected all tolerances to be checked")
        self.assertEqual(result.keys_not_in_tol, [], "Expected all key in data to also be in tolerance")

    def test_when_list_element_does_not_match_ref(self):

        # One list value differs to reference
        data = {'a': 1.0, 'b': [1.0, 2.0]}
        ref = {'a': 1.0, 'b': [1.1, 2.0]}
        tolerance = {'a': 1.e-10, 'b': 1.e-10}

        result = ErrorFinder(data, ref, tolerance)
        error = result.get_error_value('b')
        error = np.round(error, 14)

        # Assert the key of the error
        self.assertEqual(result.get_error_keys(), ['b'], msg="Only expect error in key b")
        self.assertEqual(result.get_error_keys(last=True), ['b'], msg="No nesting of keys, so last=True returns same result")

        self.assertListEqual(list(error), [0.1], msg='Only the error in the list gets returned')
        self.assertEqual(result.unused_tol_keys, [], "Expect all tolerances to be checked")
        self.assertEqual(result.keys_not_in_tol, [], "Expect all key in data to also be in tolerance")

    def test_when_2d_array_element_does_not_match(self):

        data = {1: {'energies': np.array([[0, 0, 0], [1, 1, 1]])}, 2: {'energies': np.array([[0, 0, 0], [1, 1, 1]])}}
        ref = {1: {'energies': np.array([[0, 0, 0], [1, 1, 1]])}, 2: {'energies': np.array([[0, 0, 0], [1, 2, 1]])}}
        tolerance = {'energies': 0}

        result = ErrorFinder(data, ref, tolerance)
        entry = result.get_entry('2%energies')
        error = entry.diff
        error = np.round(error, 14)

        self.assertEqual(result.get_error_keys(), ['2%energies'], msg="Only expect error in key %2%energies")
        self.assertEqual(result.get_error_keys(last=True), ['energies'], msg="Most-nested key is 'energies'")

        self.assertListEqual(list(error), [1], msg='Only the error in the list gets returned')
        self.assertEqual(result.unused_tol_keys, [], "Expected all tolerances to be checked")
        self.assertEqual(result.keys_not_in_tol, [], "Expect all key in data to also be in tolerance")

    def test_nested_key_appears_twice(self):
        """
        Test that the correct key and value is reported
        """

        data = {1: {'k_point': [0, 0, 0]}, 2: {'k_point': [1, 1, 1]}}
        ref = {1: {'k_point': [0, 1, 0]}, 2: {'k_point': [1, 1, 1]}}
        tolerance = {'k_point': 0}

        result = ErrorFinder(data, ref, tolerance)
        error = result.get_error_value('1%k_point')
        error = np.round(error, 14)

        # Assert the key of the error
        self.assertEqual(result.get_error_keys(), ['1%k_point'], msg="Only expect error in key 1%k_point")

        self.assertListEqual(list(error), [1], msg='Only the error in the list gets returned')
        self.assertEqual(result.unused_tol_keys, [], "Expected all tolerances to be checked")
        self.assertEqual(result.keys_not_in_tol, [], "Expect all key in data to also be in tolerance")

    def test_when_nested_dict_values_do_not_match_ref(self):

        data = {'a': 1.0, 'b': {'c': 1, 'd': 1.1, 'e': 'hello'}}
        ref = {'a': 1.0, 'b': {'c': 2, 'd': 5.1, 'e': 'goodbye'}}
        tolerance = {'a': 1.e-10, 'c': 0, 'd': 1.e-10, 'e': ''}

        result = ErrorFinder(data, ref, tolerance)

        self.assertEqual(len(result.errors), 3, msg='Three errors in data')
        self.assertEqual(result.get_error_keys(), ['b%c', 'b%d', 'b%e'], "Concatenated tolerance keys with special character")

        error_1 = result.get_error_value('b%c')
        error_2 = result.get_error_value('b%d')
        error_3 = result.get_error_value('b%e')

        error_1 = np.round(error_1, 14)
        error_2 = np.round(error_2, 14)

        self.assertEqual(error_1, 1., "Error in c")
        self.assertEqual(error_2, 4., "Error in d")
        self.assertEqual(error_3, 'strings differ', "Error in e")

        self.assertEqual(result.unused_tol_keys, [], "Expected all tolerances to be checked")
        self.assertEqual(result.keys_not_in_tol, [], "Expect all key in data to also be in tolerance")

    def test_when_many_values_do_not_match_ref(self):

        # Nested list and dict
        data = {'a': {'b': 1, 'c': 1.1, 'd': 'hello', 'e': {'f': 5.5}, 'g': [7, 8.9, 1.e-1]}}
        ref = {'a': {'b': 1, 'c': 1.1, 'd': 'hello', 'e': {'f': 8.3}, 'g': [8, 9.0, 1.e-2]}}
        tolerance = {'b': 0, 'c': 1.e-10, 'd': '', 'f': 1.e-10, 'g': 1.e-6}

        result = ErrorFinder(data, ref, tolerance)

        self.assertEqual(result.get_error_keys(), ['a%e%f', 'a%g'])
        error_1 = result.get_error_value('a%e%f')
        error_2 = result.get_error_value('a%g')

        error_1 = np.round(error_1, 14)
        error_2 = np.round(error_2, 14)

        self.assertEqual(error_1, 2.8, "Error in f")
        self.assertListEqual(list(error_2), [1., 0.1, 9.e-2], "Error in g")

        self.assertEqual(result.unused_tol_keys, [], "Expected all tolerances to be checked")
        self.assertEqual(result.keys_not_in_tol, [], "Expect all key in data to also be in tolerance")
