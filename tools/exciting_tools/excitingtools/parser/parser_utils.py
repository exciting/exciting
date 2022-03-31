"""Utilities for parsers
"""
from typing import Callable
import os
import xml.etree.ElementTree as ET


def generic_parser(file_name: str, parser_func: Callable[[str], dict]) -> dict:
    """
    Generic parser provides a wrapper for file IO.

    Note, this could be refactored as a decorator.

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


def xml_root(func: Callable):
    """ Decorate XML parsers, enabling the developer to pass
    an XML file name, XML string or ElementTree.Element as input
    and return the XML root.
    """
    def modified_func(input: str):
        # Element
        if isinstance(input, ET.Element):
            return func(input)

        # File name
        try:
            tree = ET.parse(input)
            root = tree.getroot()
            return func(root)
        except (FileNotFoundError, OSError):
            pass

        # XML string
        try:
            root = ET.fromstring(input)
            return func(root)
        except ET.ParseError:
            raise ValueError(f'Input string neither an XML file, '
                             f'nor valid XML: {input}')

    return modified_func
