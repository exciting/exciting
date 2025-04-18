"""Classes and functions to aid unit testing."""

import builtins
import pathlib
from types import ModuleType

real_import = builtins.__import__


class MockFile:
    """Single class for testing parsers that require either:
       * File.
       * String contents of file.

    Usage:
    @pytest.fixture
    def file_mock(tmp_path):
       file = tmp_path / _file_name
       file.write_text(string_contents)
       return MockFile(file, string_contents)
    """

    def __init__(self, file: pathlib.Path, string: str):
        # File object
        self.file = file
        # File contents
        self.string = string
        # Name prepended by path
        self.full_path = self.file.as_posix()


def import_error_monty(name: str, *args, **kwargs) -> ModuleType:
    """Mock that monty is not available.

    :param name: module to import.
    :param args: import args
    :param kwargs: import kwargs
    :return: the imported module
    """
    if name.startswith("monty"):
        raise ModuleNotFoundError("Simulating missing monty")
    return real_import(name, *args, **kwargs)
