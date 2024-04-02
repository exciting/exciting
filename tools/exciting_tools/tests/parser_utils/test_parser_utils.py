"""Tests for parser utils."""

import pytest

from excitingtools.parser_utils.parser_utils import (
    convert_single_entry,
    convert_string_dict,
    json_convert,
    standardise_fortran_exponent,
)


@pytest.mark.parametrize(
    "test_input,expected",
    [
        ("true", True),
        ("false", False),
        ("fromscratch", "fromscratch"),
        ("skip", "skip"),
        ("3", 3),
        ("-1", -1),
        ("2.34", 2.34),
        ("5.1e3", 5100),
    ],
)
def test_convert_json(test_input, expected):
    assert json_convert(test_input) == expected


@pytest.mark.parametrize("test_input,expected", [("23.2D-1", 2.32), ("1q0", 1)])
def test_convert_fortran_exponent(test_input, expected):
    assert standardise_fortran_exponent(test_input, return_as_str=False) == expected


@pytest.mark.parametrize(
    "test_input,expected",
    [
        ("true", True),
        ("skip", "skip"),
        ("3", 3),
        ("2.34", 2.34),
        ("true false false", [True, False, False]),
        ("1.3 2.3e4 3 90", [1.3, 2.3e4, 3, 90]),
        ("4 4 3", [4, 4, 3]),
        ("1 2", [1, 2]),
        ("ab cd", "ab cd"),
    ],
)
def test_convert_single_entry(test_input, expected):
    assert convert_single_entry(test_input) == expected


def test_convert_single_entry_to_int():
    converted_int = convert_single_entry("3")
    assert converted_int == 3
    # need to test really the type to detect the error
    assert isinstance(converted_int, int)


def test_convert_string_dict():
    string_dict = {"a": "1", "b": "b", "c": "false", "d": "3 3 3"}
    ref_dict = {"a": 1, "b": "b", "c": False, "d": [3, 3, 3]}
    assert convert_string_dict(string_dict) == ref_dict
