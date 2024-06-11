"""General utils for exciting scripts."""

import numpy as np
from typing import Tuple, List, TypeVar, Union

# Types for static type checking to support maintaining type consistency. For example, the function
# `sort_lists_by_first_list` should return a tuple of lists with elements belonging to the same type as the elements in
# the lists passed as arguments.
T1 = TypeVar("T1")
T2 = TypeVar("T2")


def sort_lists_by_first_list(first_list: List[T1], second_list: List[T2]) -> Tuple[List[T1], List[T2]]:
    """ Sorts two lists, using the first list as reference

    :param first_list: first list to be sorted, used as reference
    :param second_list: second list to be sorted, uses first list as reference
    :return: sorted lists
    """
    first_len = len(first_list)
    second_len = len(second_list)
    assert first_len == second_len, f"Both lists should have the same length, not {first_len} and {second_len}."

    sorted_indices = np.argsort(first_list)
    sorted_first_list = [first_list[x] for x in sorted_indices]
    sorted_second_list = [second_list[x] for x in sorted_indices]

    return sorted_first_list, sorted_second_list

def get_decimal_decomposition(number: float) -> Tuple[float, int]:
    """Decompose the number into mantissa and exponent.

    :param number: input number
    :return: tuple with shifted number (only one leading digit before the decimal point) and exponent
    """
    log_num = np.log10(abs(number))
    exponent = int(log_num)
    shifted_number = 10 ** (log_num - exponent)
    while shifted_number < 1:
        shifted_number *= 10
        exponent -= 1
    return shifted_number, exponent


def get_prettified_scientific_notation(number: float, unit: Union[str, None] = None) -> str:
    """Decompose the number into mantissa and exponent and produce formatted string.

    :param number: input number
    :param unit: unit of the number
    :return: prettified string representation
    """
    shifted_number, exponent = get_decimal_decomposition(number)

    exponent_string = f"$10^{{{exponent}}}$"
    if exponent == 0:
        exponent_string = ""
    elif exponent == 1:
        exponent_string = r"$10\,$"

    sign = "+" if number >= 0 else "\u2013"
    representation = rf"${sign}{shifted_number:6.4f}\cdot${exponent_string}"
    if unit is None:
        return representation
    return representation + f"[{unit}]"
