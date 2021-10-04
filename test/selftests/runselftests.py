"""
Unit tests for the test suite. Run the tests:
cd test/
python3 selftests/runselftests.py
"""
import unittest
import json
import os

from modules.tester.report import Test
from modules.tolerance.modify_tolerance import extend_tolerance_file_keys, remove_tolerance_file_keys


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
        test = Test(name='NoFailures', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 0, 'format': 0, 'float': 0, 'array': 0, 'string': 0, 'integer': 0, 'reference': 0}))

    def test_compare_data_FormatFails(self):
        data = {'a': 1.0, 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        reference = {'b': 1.0, 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        test = Test(name='FormatFail. Difference in first key', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 1, 'format': 1, 'float': 0, 'array': 0, 'string': 0, 'integer': 0, 'reference': 0}))

        data = {'a': 1.0, 'b': [1.0, 2.0, 4.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        reference = {'a': 1.0, 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        test = Test(name='FormatFail. Difference in b value', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 1, 'format': 1, 'float': 0, 'array': 0, 'string': 0, 'integer': 0, 'reference': 0}))

        data = {'a': 1.0, 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], "a"], 'f': ['1', '2']}
        reference = {'a': 1.0, 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        test = Test(name='FormatFail. Difference in e value', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 1, 'format': 1, 'float': 0, 'array': 0, 'string': 0, 'integer': 0, 'reference': 0}))

        data = {'a': 'a', 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        reference = {'a': 1.0, 'b': [1.0, 2, 0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        test = Test(name='FormatFail. Difference in a value', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 1, 'format': 1, 'float': 0, 'array': 0, 'string': 0, 'integer': 0, 'reference': 0}))

    def test_compare_data_ListFails(self):
        data = {'a': 5.0, 'b': [1.0, 10.0], 'c': {'a1': 1.0}, 'd': "blu", 'e': [[1, 2], 3], 'f': ['1', '2']}
        reference = {'a': 5.0, 'b': [1.0, 2.0], 'c': {'a1': 1.0}, 'd': "blu", 'e': [[1, 2], 3], 'f': ['1', '2']}
        test = Test(name='ListFail', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 1, 'format': 0, 'float': 0, 'array': 1, 'string': 0, 'integer': 0, 'reference': 0}))

    def test_compare_data_FloatFails(self):
        data = {'a': 5.0, 'b': [1.0, 10.0], 'c': {'a1': 1.0}, 'd': "blu", 'e': [[1, 2], 3], 'f': ['1', '2']}
        reference = {'a': 1.0, 'b': [1.0, 10.0], 'c': {'a1': 1.0}, 'd': "blu", 'e': [[1, 2], 3], 'f': ['1', '2']}
        test = Test(name='FloatFail', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 1, 'format': 0, 'float': 1, 'array': 0, 'string': 0, 'integer': 0, 'reference': 0}))

    def test_compare_data_StringFails(self):
        data = {'a': 5.0, 'b': [1.0, 10.0], 'c': {'a1': 1.0}, 'd': "blu", 'e': [[1, 2], 3], 'f': ['1', '2']}
        reference = {'a': 5.0, 'b': [1.0, 10.0], 'c': {'a1': 1.0}, 'd': "bla", 'e': [[1, 2], 3], 'f': ['1', '2']}
        test = Test(name='StringFail', tolFloat=1.0, tolMSE=1.0, tolValuePair=self.tol_value_pair)
        test.evaluate(data, reference, self.test_dir)
        self.assertEqual(test.count_failures(),
                         dict({'total': 1, 'format': 0, 'float': 0, 'array': 0, 'string': 1, 'integer': 0, 'reference': 0}))


class TestModifyTolerances(unittest.TestCase):
    """
    Test routines extend_tolerance_file_keys and remove_tolerance_file_keys
    """

    initial_tolerances = {"info_out":
                              {'Lattice vectors (cartesian)': 1.e-8,
                               'Reciprocal lattice vectors (cartesian)': 1.e-8,
                               'Unit cell volume': 1.e-6,
                               'Number of crystal symmetries': 0
                               }
                          }

    file_name = 'gs_tols.json'

    def setUp(self) -> None:
        """
        Create reference file for tests to read
        """
        with open(self.file_name, "w") as fid:
            json.dump(self.initial_tolerances, fid, indent=2)
        return

    def test_extend_tolerance_file_keys(self):
        """
        Interacts with the file system (writes, reads, then deletes)

        Default behaviour tested, hence no lambda passed as an optional arg
        """

        new_entries = {"info_out":
                           {'some_new_key': 1.e-5}
                       }

        final_tolerances = {"info_out":
                                {'Lattice vectors (cartesian)': 1.e-8,
                                 'Reciprocal lattice vectors (cartesian)': 1.e-8,
                                 'Unit cell volume': 1.e-6,
                                 'Number of crystal symmetries': 0,
                                 'some_new_key': 1.e-5
                                 }
                            }

        # Modify tolerance file
        missed_files = extend_tolerance_file_keys([self.file_name], new_entries)
        self.assertTrue(not missed_files, "Expect " + self.file_name + " to be found and modified")

        # Check file now contains new_entries
        with open(self.file_name) as fid:
            actual_tolerances = json.load(fid)

        self.assertEqual(actual_tolerances, final_tolerances,
                         "Expect actual_tolerances to contain old data plus {new key:tol} entry")

        # Check original file was backed-up and is unchanged
        original_file = 'orig_' + self.file_name
        self.assertTrue(os.path.isfile(original_file), "Expect a copy of the original file to be made")

        with open(original_file) as fid:
            original_tolerances = json.load(fid)

        self.assertEqual(original_tolerances, self.initial_tolerances,
                         "Expect original file to contain initial_tolerances")

        # Clean files
        os.remove(self.file_name)
        os.remove(original_file)

    def test_extend_tolerance_file_keys_without_backup(self):
        """
        Instead of the default behaviour, overwrite the original tolerance file with the modified file
        """

        new_entries = {"info_out":
                           {'some_new_key': 1.e-5}
                       }

        final_tolerances = {"info_out":
                                {'Lattice vectors (cartesian)': 1.e-8,
                                 'Reciprocal lattice vectors (cartesian)': 1.e-8,
                                 'Unit cell volume': 1.e-6,
                                 'Number of crystal symmetries': 0,
                                 'some_new_key': 1.e-5
                                 }
                            }

        # Extend tolerance file
        missed_files = extend_tolerance_file_keys([self.file_name], new_entries, rename_original=lambda name: name)
        self.assertTrue(not missed_files, "Expect " + self.file_name + " to be found and modified")

        # Check file now contains new_entries
        with open(self.file_name) as fid:
            actual_tolerances = json.load(fid)

        self.assertEqual(actual_tolerances, final_tolerances,
                         "Expect actual_tolerances to contain old data plus {new key:tol} entry")

        # Check that a back-up was not created
        original_file = 'orig_' + self.file_name
        self.assertTrue(not os.path.isfile(original_file), "Do not expect a back-up of the original file to be made")

        # Clean file
        os.remove(self.file_name)

    def test_remove_tolerance_file_keys(self):
        """
        Interacts with the file system (writes, reads, then deletes)

        Default behaviour tested, hence no lambda passed as an optional arg
        """

        entry_to_remove = {"info_out":
                               ['Number of crystal symmetries']
                           }

        final_tolerances = {"info_out":
                                {'Lattice vectors (cartesian)': 1.e-8,
                                 'Reciprocal lattice vectors (cartesian)': 1.e-8,
                                 'Unit cell volume': 1.e-6
                                 }
                            }

        missed_files = remove_tolerance_file_keys([self.file_name], entry_to_remove)
        self.assertTrue(not missed_files, "Expect " + self.file_name + " to be found and modified")

        # Check file contains original tolerances minus entry_to_remove
        with open(self.file_name) as fid:
            actual_tolerances = json.load(fid)

        self.assertEqual(actual_tolerances, final_tolerances,
                         "Expect actual_tolerances to contain old data minus {Number of crystal symmetries:tol}")

        # Check original file was backed-up and is unchanged
        original_file = 'orig_' + self.file_name
        self.assertTrue(os.path.isfile(original_file), "Expect a copy of the original file to be made")

        with open(original_file) as fid:
            original_tolerances = json.load(fid)

        self.assertEqual(original_tolerances, self.initial_tolerances,
                         "Expect original file to contain initial_tolerances")

        # Clean files
        os.remove(self.file_name)
        os.remove(original_file)

    def test_remove_tolerance_file_keys_without_backup(self):
        """
        Interacts with the file system (writes, reads, then deletes)

        Instead of the default behaviour, overwrite the original tolerance file with the modified file
        """

        entry_to_remove = {"info_out":
                               ['Number of crystal symmetries']
                           }

        final_tolerances = {"info_out":
                                {'Lattice vectors (cartesian)': 1.e-8,
                                 'Reciprocal lattice vectors (cartesian)': 1.e-8,
                                 'Unit cell volume': 1.e-6
                                 }
                            }

        # Remove entry in tolerance file
        missed_files = remove_tolerance_file_keys([self.file_name], entry_to_remove, rename_original=lambda name: name)
        self.assertTrue(not missed_files, "Expect " + self.file_name + " to be found and modified")

        # Check file contains original tolerances minus entry_to_remove
        with open(self.file_name) as fid:
            actual_tolerances = json.load(fid)

        self.assertEqual(actual_tolerances, final_tolerances,
                         "Expect actual_tolerances to contain old data plus {new key:tol} entry")

        # Check that a back-up was not created
        original_file = 'orig_' + self.file_name
        self.assertTrue(not os.path.isfile(original_file), "Do not expect a back-up of the original file to be made")

        # Clean file
        os.remove(self.file_name)


if __name__ == '__main__':
    unittest.main()
