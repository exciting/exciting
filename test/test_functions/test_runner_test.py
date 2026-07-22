"""Tests for helper logic in runner/test.py."""

from ..src.runner.test import get_test_main_output_file_name


def test_get_test_main_output_file_name_uses_resolved_properties_output():
    output_files = ["EPSILON_11.OUT", "EPSILON_12.OUT", "CHI_111.OUT"]
    assert get_test_main_output_file_name("properties", output_files) == "EPSILON_11.OUT"


def test_get_test_main_output_file_name_keeps_info_out_when_present():
    output_files = ["INFO.OUT", "MBD.OUT"]
    assert get_test_main_output_file_name("properties", output_files) == "INFO.OUT"


def test_get_test_main_output_file_name_non_properties_uses_default():
    output_files = ["evalcore.xml", "eigval.xml"]
    assert get_test_main_output_file_name("groundstate", output_files) == "INFO.OUT"
