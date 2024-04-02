"""Decorators and wrappers for parser functions."""

import pathlib
import xml.etree.ElementTree as ET
from typing import Callable, Optional, Union

from excitingtools.utils.dict_utils import __container_converter


def return_file_string(file_name: Union[str, pathlib.Path]) -> str:
    """Given a file name, return the file contents as a string.

    :param file_name: File name.
    :return file_string: File contents string.
    """
    file_name_ = file_name

    if isinstance(file_name_, str):
        file_name_ = pathlib.Path(file_name_)

    if not file_name_.exists():
        raise FileNotFoundError(f"{file_name_} not found")

    return file_name_.read_text()


def file_handler(file_name: Union[str, pathlib.Path], parser_func: Callable[[str], dict]) -> dict:
    """Provide a wrapper for file IO.

    :param file_name: File name or Path object
    :param Callable[[str], dict] parser_func: Parser function, which expects a parsed
     string as its only input and returns a dictionary.
    :return: dict data: Dictionary of parsed data, with values converted from strings.
    """
    file_string = return_file_string(file_name)
    return parser_func(file_string)


def accept_file_name(parser: Callable):
    """Decorate parsers that accept string contents, such that they take file names instead."""

    def modified_func(file_name: Union[str, pathlib.Path]):
        """Wrapper.
        param: file_name: File name.
        """
        file_string = return_file_string(file_name)
        return parser(file_string)

    return modified_func


def set_return_values(parser: Callable[[str], dict]) -> Callable[[str], dict]:
    """Mutate the values of a parsed dictionary to return
    appropriate types, rather than strings.
    """

    def modified_exciting_parser(full_file_name: str) -> dict:
        """Wrapper.
        :param full_file_name: File name.
        :return: converted data
        """
        data = parser(full_file_name)
        return __container_converter(data)

    return modified_exciting_parser


def xml_root(func: Callable):
    """Decorate XML parsers, enabling the developer to pass
    an XML file name, XML string or ElementTree.Element as input
    and return the XML root.
    """
    function_selection = {type(None): lambda x, _: func(x), str: lambda x, y: func(x, y)}

    def modified_func(input: str, tag: Optional[str] = None):
        # Element
        if isinstance(input, ET.Element):
            return function_selection[type(tag)](input, tag)

        # File name
        try:
            tree = ET.parse(input)
            root = tree.getroot()
            return function_selection[type(tag)](root, tag)
        except (FileNotFoundError, OSError):
            pass

        # XML string
        try:
            root = ET.fromstring(input)
            return function_selection[type(tag)](root, tag)
        except (ET.ParseError, TypeError):
            raise ValueError(f"Input string neither an XML file, nor valid XML: {input}")

    return modified_func
