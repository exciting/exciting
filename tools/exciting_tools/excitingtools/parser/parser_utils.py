"""
Utilities for parsers
"""
from typing import Callable
import os


# TODO(Alex) This should be a decorator that takes an argument
def generic_parser(file_name: str, parser_func: Callable[[str], dict]) -> dict:
    """
    Generic parser provides a wrapper for file IO.

    :param str file_name: Name of file to open and parse.
    :param Callable[[str], dict] parser_func: Parser function, which expects a parsed
     string as its only input and returns a dictionary.

    :return: dict data: Dictionary of parsed data, with values converted from strings.
    """
    if not os.path.exists(file_name):
        raise OSError('File path not valid:', file_name)

    with open(file_name) as f:
        file_string = f.read()

    return parser_func(file_string)
