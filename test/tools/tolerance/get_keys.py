"""
Module to print the unique keys found in the template tolerance files.

Given a results dictionary, as returned by an exciting parser, this returns the lowest-level set of nested keys.
These are the keys with corresponding values, that are compared by the test suite.
As such, each of these keys should have an associated tolerance.

To print the keys associated with a parsed exciting output file, type:

```
  python3 get_keys.py -file path/2/file
```

"""
import argparse

from excitingtools.parser import parser_chooser
from excitingtools.dict_utils import get_hashable_entries


def parse_input_args() -> dict:
    """
    Parse command line inputs

    :return dict: Dictionary of parsed inputs
    """
    parser = argparse.ArgumentParser(description="Module to generate keys consistent with parsing exciting output files")

    parser.add_argument('-file',
                        help='Given an exciting output file, get and print the keys associated with '
                             'parsing it (parser-specific)',
                        dest='file_names',
                        type=str,
                        nargs='+'
                        )

    return parser.parse_args().__dict__


def print_tolerance_keys(file: str):
    """
    Given some exciting output file, print all unique keys.

    For use in generating the keys used in each respective python file
    containing default tolerances.

    :param file: file name, prepended by the file path
    """
    output = parser_chooser(file)

    unique_keys = set([key for key, value in get_hashable_entries(output)])

    for unique_key in unique_keys:
        print(" " + unique_key)


if __name__ == "__main__":
    print("Printing keys associated with parsed exciting output files.")

    inputs = parse_input_args()

    for name in inputs['file_names']:
        print("File:", name.split('/')[-1])
        print_tolerance_keys(name + '\n')
        print("\n\n")
