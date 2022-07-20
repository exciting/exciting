"""General utility functions. Typically conversion/type-checking.
"""
from typing import Union, List, Optional, Callable
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


def convert_to_literal(input_str: str) -> Union[int, float, None]:
    """
    If possible, convert string to an int or float

    For example:
      convert_to_literal('1.1') returns 1.1
      convert_to_literal('1') returns 1
      convert_to_literal('1.0') returns 1.0

    :param str input_str: Input string
    :return Union[int, float] x: Numerical literal of x, else None.
    """
    try:
        integer_string = int(input_str)
        return integer_string
    except ValueError:
        try:
            float_string = float(input_str)
            return float_string
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


def list_to_str(mylist: list, modifier: Optional[Callable] = None) -> str:
    """ Convert a list to a string
    """
    if modifier is None:
        modifier = lambda x: x
    return "".join(modifier(str(xyz)) + ' ' for xyz in mylist).strip()


def string_to_bool(string: str) -> bool:
    """ Convert string representation of true/false to True/False.

    :param string: String
    :return bool
    """
    if string.lower() == 'true':
        return True
    elif string.lower() == 'false':
        return False
    else:
        raise ValueError()



