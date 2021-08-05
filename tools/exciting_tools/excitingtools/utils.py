"""
General utility functions
"""
from typing import Union

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
