"""
Test functions in set_inputs.py
"""
import pathlib

from ..src.runner.set_inputs import input_files_for_tests


def test_input_files_for_tests(tmp_path):
    """
    Test if input files can be found from inspecting the file system.

    Note, `test_case_species_files` is the only routine that actually interacts with
    the file system.
    """
    test_name: pathlib.Path = tmp_path / "test_Ar"

    # Mock directories based on test name
    test_name.mkdir()
    ref_directory = test_name / "ref"
    ref_directory.mkdir()
    assert ref_directory.is_dir(), "tmp_path directory does not exist"

    # Mock input files
    input_xml = ref_directory / "input.xml"
    species_xml = ref_directory / "Ar.xml"
    input_xml.write_text('Arbitrary data')
    species_xml.write_text('Arbitrary data')

    test_name_str = test_name.as_posix()
    input_files = input_files_for_tests([test_name_str], subdirectory='ref')

    assert list(input_files.keys()) == [test_name_str], "Input file keys should equal the test names"
    assert input_files[test_name_str] == ['Ar.xml', 'input.xml'], "Found erroneous input files"
