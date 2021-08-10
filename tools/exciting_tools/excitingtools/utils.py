"""
General utility functions
"""
from typing import Union, List
import re

def can_be_float(value) -> bool:
    """
    Check if a value can be interpreted as a float

    :param value: Input
    :return bool: Value can be interpreted as a float
    """
    try:
        float(value)
        return True
    except ValueError:
        return False


def convert_to_literal(s: str) -> Union[int, float]:
    """
    If possible, convert string to an int or float

    For example:
      convert_to_literal('1.1') returns 1.1
      convert_to_literal('1') returns 1
      convert_to_literal('1.0') returns 1.0

    :param str s: Input string
    :return Union[int, float] x: Numerical literal of x, else None.
    """
    try:
        x = int(s)
        return x
    except ValueError:
        try:
            x = float(s)
            return x
        except ValueError:
            return None


def get_new_line_indices(string: str) -> List[int]:
    """
    Given a string, return the indices that correspond to the
    start of new lines.

    For example,
     line = get_new_line_indices(string)
     # First line
     string[line[0]: line[1]]
     # 6th line
     string[line[5]: line[6]]

    :param str string: Input string
    :return List[int] indices: List of indices corresponding to
     new lines in string.
    """
    indices = [0]
    indices += [m.start() + 1 for m in re.finditer('\n', string)]
    return indices
