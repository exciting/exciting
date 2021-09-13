import unittest
import numpy as np

from tools.tester.compare import hashable_list, compare_reference_with_target


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


class TestCompareReferenceWithTarget(unittest.TestCase):
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
        errors, unused_tol_keys = compare_reference_with_target(data, data, tolerance)
        self.assertEqual(errors, {}, msg='Each hashable data type in compare_reference_with_target')
        self.assertEqual(unused_tol_keys, [], "Expected all tolerances to be checked")

        # Value = string, with lr whitespace ignored
        data = {'a': '  hello  '}
        ref = {'a': 'hello'}
        errors, unused_tol_keys = compare_reference_with_target(data, ref, {'a': ''})
        self.assertEqual(errors, {}, msg="Left and right whitespace ignored in string comparisons")
        self.assertEqual(unused_tol_keys, [], "Expected all tolerances to be checked")

        # Value = list
        data = {'a': 1.0, 'b': [1.0, 2.0]}
        tolerance = {'a': 1.e-10, 'b': 1.e-10}
        errors, unused_tol_keys = compare_reference_with_target(data, data, tolerance)
        self.assertEqual(errors, {}, msg="Value in target dictionary = list")
        self.assertEqual(unused_tol_keys, [], "Expected all tolerances to be checked")

        # Value = dict
        data = {'a': 1.0, 'b': {'c': 1, 'd': 1.1, 'e': 'hello'}}
        tolerance = {'a': 1.e-10, 'c': 0, 'd': 1.e-10, 'e': ''}
        errors, unused_tol_keys = compare_reference_with_target(data, data, tolerance)
        self.assertEqual(errors, {}, msg="Value in target dictionary = dictionary")
        self.assertEqual(unused_tol_keys, [], "Expected all tolerances to be checked")

        # Value = Nested list and dict
        data = {'a': {'b': 1, 'c': 1.1, 'd': 'hello', 'e': {'f': 5.5}, 'g': [7, 8.9, 1.e-9]}}
        tolerance = {'b': 0, 'c': 1.e-10, 'd': '', 'f': 1.e-10, 'g': 1.e-10}
        errors, unused_tol_keys = compare_reference_with_target(data, data, tolerance)
        self.assertEqual(errors, {}, msg="Values in target dictionary contain a mix of dict and list")
        self.assertEqual(unused_tol_keys, [], "Expected all tolerances to be checked")

    def test_ignored_tolerance_keys(self):
        data = {'a': 1, 'b': 1.1}
        tolerance = {'a': 0, 'b': 1.e-10, 'c': np.infty}
        errors, unused_tol_keys = compare_reference_with_target(data, data, tolerance)
        self.assertEqual(unused_tol_keys, ['c'], msg="Keys of unevaluated tolerances")
        self.assertEqual(errors, {}, msg="Expect errors to be empty")

    def test_string_tolerances(self):
        data = {'a': 'Some string'}
        tolerance = {'a': 'arbitrary string tolerance'}
        errors, unused_tol_keys = compare_reference_with_target(data, data, tolerance)
        self.assertEqual(unused_tol_keys, [], msg="Expected key with string tolerance to be accessed")
        self.assertEqual(errors, {}, msg="String tolerance values are not evaluated. Only the keys are used")

    def test_when_values_do_not_match_ref(self):
        """
        Test when the test/target values do not match the reference, to within tolerance.
        """

        # Each hashable data type
        data = {'a': 1, 'b': 1.1, 'c': 'hello'}
        ref = {'a': 2, 'b': 1.2, 'c': 'goodbye'}
        tolerance = {'a': 0, 'b': 1.e-10, 'c': ''}

        errors, unused_tol_keys = compare_reference_with_target(data, ref, tolerance)
        errors = round_values(errors, dp=14)

        self.assertEqual(errors, {'a': 1, 'b': 0.1, 'c': 'strings differ'}, msg="Expect errors in keys a, b and c")
        self.assertEqual(unused_tol_keys, [], "Expected all tolerances to be checked")

        # All list values in data differ to reference
        data = {'a': 1.0, 'b': [1.0, 2.0]}
        ref = {'a': 1.0, 'b': [1.1, 2.5]}
        tolerance = {'a': 1.e-10, 'b': 1.e-10}

        errors, unused_tol_keys = compare_reference_with_target(data, ref, tolerance)
        errors = round_values(errors, dp=14)

        self.assertEqual(list(errors.keys()), ['b'], msg="Only expect error in key b")
        self.assertListEqual(list(errors['b']), [0.1, 0.5], msg="All values in list differ to ref")
        self.assertEqual(unused_tol_keys, [], "Expected all tolerances to be checked")

        # One list value differs to reference
        data = {'a': 1.0, 'b': [1.0, 2.0]}
        ref = {'a': 1.0, 'b': [1.1, 2.0]}
        tolerance = {'a': 1.e-10, 'b': 1.e-10}

        errors, unused_tol_keys = compare_reference_with_target(data, ref, tolerance)
        errors = round_values(errors, dp=14)

        self.assertEqual(list(errors.keys()), ['b'])
        self.assertListEqual(list(errors['b']), [0.1, 0.0], msg='Only one error in the list but all diffs get returned')
        self.assertEqual(unused_tol_keys, [], "Expected all tolerances to be checked")

        # Some value = dict
        data = {'a': 1.0, 'b': {'c': 1, 'd': 1.1, 'e': 'hello'}}
        ref = {'a': 1.0, 'b': {'c': 2, 'd': 5.1, 'e': 'goodbye'}}
        tolerance = {'a': 1.e-10, 'c': 0, 'd': 1.e-10, 'e': ''}

        errors, unused_tol_keys = compare_reference_with_target(data, ref, tolerance)
        errors = round_values(errors, dp=14)

        self.assertEqual(errors, {'c': 1, 'd': 4.0, 'e': 'strings differ'}, msg='Expect values of c, d and e to differ')
        self.assertEqual(unused_tol_keys, [], "Expected all tolerances to be checked")

        # Nested list and dict
        data = {'a': {'b': 1, 'c': 1.1, 'd': 'hello', 'e': {'f': 5.5}, 'g': [7, 8.9, 1.e-1]}}
        ref = {'a': {'b': 1, 'c': 1.1, 'd': 'hello', 'e': {'f': 8.3}, 'g': [8, 9.0, 1.e-2]}}
        tolerance = {'b': 0, 'c': 1.e-10, 'd': '', 'f': 1.e-10, 'g': 1.e-6}

        errors, unused_tol_keys = compare_reference_with_target(data, ref, tolerance)
        errors = round_values(errors, dp=14)

        self.assertEqual(list(errors.keys()), ['f', 'g'])
        self.assertEqual(errors['f'], 2.8)
        self.assertListEqual(list(errors['g']), [1, 0.1,  0.09])
        self.assertEqual(unused_tol_keys, [], "Expected all tolerances to be checked")
