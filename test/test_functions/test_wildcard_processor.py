"""
Test the function wildcard_processor
To run with pytest:
```
cd test
pytest -s selftests/test_wildcard_processor.py
```
"""
import unittest
from ..src.utilities.wildcard_processor import wildcard_processor


class test_wild_card_processor(unittest.TestCase):
    """
    Test wildcard_processor.
    """

    def test_no_wildcard(self):
        """
        Test for a string without wildcard that only start of string and end of string
        symbols are added, and that '.' is replaced by '[.]'.
        """
        # String without '.'
        string = r'ased_asdneudj'
        test_string = wildcard_processor(string)
        ref_string = r'^ased_asdneudj$'
        self.assertEqual(test_string, ref_string, 'Expect that only start and end symbols are included.')

        # String with '.'
        string = r'ased.asdneudj'
        test_string = wildcard_processor(string)
        ref_string = r'^ased[.]asdneudj$'
        self.assertEqual(test_string, ref_string,
                         'Expect that only start and end symbols are included, and "." is replaced by "[.]".')

    def test_letter_wildcard(self):
        """
        Test for a string with letter wildcard '&'.
        """
        string = r'ased_&&_asd&udj'
        test_string = wildcard_processor(string)
        ref_string = r'^ased_[a-zA-Z][a-zA-Z]_asd[a-zA-Z]udj$'
        self.assertEqual(test_string, ref_string, 'Expect that each "&" is replaced with letter regex.')

    def test_digit_wildcard(self):
        """
        Test for a string with digit wildcard '#'.
        """
        string = r'ased_##_asd#udj'
        test_string = wildcard_processor(string)
        ref_string = r'^ased_[\d][\d]_asd[\d]udj$'
        self.assertEqual(test_string, ref_string, 'Expect that each "#" is replaced with digit regex".')

    def test_character_wildcard(self):
        """
        Test for a string with character wildcard '*'.
        """
        string = r'ased_**_asd*udj'
        test_string = wildcard_processor(string)
        ref_string = r'^ased_[a-zA-Z\d][a-zA-Z\d]_asd[a-zA-Z\d]udj$'
        self.assertEqual(test_string, ref_string, 'Expect that each "&" is replaced with word character regex.')

    def test_word_wildcard(self):
        """
        Test for a string with word wildcard '??'.
        """
        string = r'ased_??_asd??udj'
        test_string = wildcard_processor(string)
        ref_string = r'^ased_[\w-]+_asd[\w-]+udj$'
        self.assertEqual(test_string, ref_string, 'Expect that each "??" is replaced with string regex.')

    def test_mixed_wildcards(self):
        """
        Test for a string with mixed wildcards.
        """
        string = r'a*ed_??_asd&u_#dj.**.??.out'
        test_string = wildcard_processor(string)
        ref_string = r'^a[a-zA-Z\d]ed_[\w-]+_asd[a-zA-Z]u_[\d]dj[.][a-zA-Z\d][a-zA-Z\d][.][\w-]+[.]out$'
        self.assertEqual(test_string, ref_string, 'Expect that each wild card is replaced with the correct regex.')
