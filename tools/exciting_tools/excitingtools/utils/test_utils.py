""" Classes and functions to aid unit testing.
"""
import pathlib


class MockFile:
    """ Single class for testing parsers that require either:
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
