"""
Unit tests for the test suite. Run the tests:
cd test/
python3 test_functions/test_old_class_test.py

TODO(A/B/H) Issue 99. Remove Functions Associated with XML Tols
"""
import unittest

# Rename on import so pytest does not attempt to capture the class under test
from ..modules.tester.report import Test as FailureLogger


class TestCompare(unittest.TestCase):
    """
    Given some reference data and some actual data, test whether the method countFailures() of Test
    returns a dictionary logging the differences between the two.
    """

    tol_value_pair = {'tol': '0', 'value': ''}
    test_dir = 'test_dir'

    def test_compare_data_NoFails(self):
        data = {'a': 1.0, 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        reference = {'a': 1.0, 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        test = FailureLogger(name='NoFailures', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 0, 'format': 0, 'float': 0, 'array': 0, 'string': 0, 'integer': 0, 'reference': 0}))

    def test_compare_data_FormatFails(self):
        data = {'a': 1.0, 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        reference = {'b': 1.0, 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        test = FailureLogger(name='FormatFail. Difference in first key', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 1, 'format': 1, 'float': 0, 'array': 0, 'string': 0, 'integer': 0, 'reference': 0}))

        data = {'a': 1.0, 'b': [1.0, 2.0, 4.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        reference = {'a': 1.0, 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        test = FailureLogger(name='FormatFail. Difference in b value', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 1, 'format': 1, 'float': 0, 'array': 0, 'string': 0, 'integer': 0, 'reference': 0}))

        data = {'a': 1.0, 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], "a"], 'f': ['1', '2']}
        reference = {'a': 1.0, 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        test = FailureLogger(name='FormatFail. Difference in e value', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 1, 'format': 1, 'float': 0, 'array': 0, 'string': 0, 'integer': 0, 'reference': 0}))

        data = {'a': 'a', 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        reference = {'a': 1.0, 'b': [1.0, 2, 0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        test = FailureLogger(name='FormatFail. Difference in a value', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 1, 'format': 1, 'float': 0, 'array': 0, 'string': 0, 'integer': 0, 'reference': 0}))

    def test_compare_data_ListFails(self):
        data = {'a': 5.0, 'b': [1.0, 10.0], 'c': {'a1': 1.0}, 'd': "blu", 'e': [[1, 2], 3], 'f': ['1', '2']}
        reference = {'a': 5.0, 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "blu", 'e': [[1, 2], 3], 'f': ['1', '2']}
        test = FailureLogger(name='ListFail', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 1, 'format': 0, 'float': 0, 'array': 1, 'string': 0, 'integer': 0, 'reference': 0}))

    def test_compare_data_FloatFails(self):
        data = {'a': 5.0, 'b': [1.0, 10.0], 'c': {'a1': 1.0}, 'd': "blu", 'e': [[1, 2], 3], 'f': ['1', '2']}
        reference = {'a': 1.0, 'b': [1.0, 10.0], 'c': {'a1': 1.0}, 'd': "blu", 'e': [[1, 2], 3], 'f': ['1', '2']}
        test = FailureLogger(name='FloatFail', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 1, 'format': 0, 'float': 1, 'array': 0, 'string': 0, 'integer': 0, 'reference': 0}))

    def test_compare_data_StringFails(self):
        data = {'a': 5.0, 'b': [1.0, 10.0], 'c': {'a1': 1.0}, 'd': "blu", 'e': [[1, 2], 3], 'f': ['1', '2']}
        reference = {'a': 5.0, 'b': [1.0, 10.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        test = FailureLogger(name='StringFail', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 1, 'format': 0, 'float': 0, 'array': 0, 'string': 1, 'integer': 0, 'reference': 0}))


if __name__ == '__main__':
    unittest.main()
