"""General parser utility functions."""

import re
from json import JSONDecodeError, loads
from typing import Any, Dict
from xml.etree import ElementTree


def find_element(root: ElementTree.Element, tag: str) -> ElementTree.Element:
    """Finds a given tag in an Element, either return the full ElementTree
    if the tag is correct or find the tag in that ElementTree.
    :param root: Element to find the tag in
    :param tag: tag to search for
    :return: the (sub)element with the correct tag
    """
    if root.tag == tag:
        return root
    return root.find(tag)


def convert_string_dict(inputs: dict) -> Dict[str, Any]:
    """Parses and converts a dictionary with string values to their actual data types.

    :param inputs: input dictionary
    :return: the converted dictionary
    """
    for key, value in inputs.items():
        inputs[key] = convert_single_entry(value)
    return inputs


def convert_single_entry(input_str: str):
    """Converts a single string. Accepts also lists separated by whitespaces.

    :param input_str: input string
    :return: the converted data
    """
    result = [json_convert(standardise_fortran_exponent(i)) for i in input_str.split()]

    if len(result) == 1:
        return result[0]
    if len(list(filter(lambda x: isinstance(x, str), result))) == len(result):
        return " ".join(result)
    return result


def json_convert(input_str: str):
    """Tries to convert a single string with no whitespaces to its actual data type
    using the json decoder (detects int, float and bool). Else returns the string as string.
    :param input_str: input string
    :return: the converted value
    """
    try:
        return loads(input_str)
    except JSONDecodeError:
        pass

    # the json standard (https://www.json.org/json-en.html) does not allow for numbers to have leading or trailing
    # decimal points, therefore we need another try to parse floats.
    try:
        return float(input_str)
    except ValueError:
        return input_str


def standardise_fortran_exponent(input_str: str, return_as_str: bool = True):
    """Tries to convert a single string representing a float value in scientific notation
    to the actual number.
    d(D) and q(Q) correspond to higher floating point precision than e(E). Replace
    them since they cannot be parsed by JSON

    :param input_str: input string
    :param return_as_str: return the number as a string in python format with e as exponential
    :return: If string can be converted to a float return it, else
    return the input string.
    """
    subbed_str = re.sub("[dDqQ]", "E", input_str, count=1)
    if subbed_str == input_str:
        return input_str

    try:
        number = float(subbed_str)
    except ValueError:
        return input_str
    if return_as_str:
        return str(number)
    return number
