"""
Wrappers for parsing with regex
"""
import re
from typing import List, Union
from math import floor

from excitingtools.utils import convert_to_literal


def parse_value_regex(file_string: str, key: str, no_colon=True) -> dict:
    """
    Match the first instance of a string (key) if present in file_string,
    and return the result in a dictionary.

    :param str file_string: Input string
    :param str key: String to match, also used as a key in the returned dictionary
    :param bool no_colon: Remove trailing colons from parsed data keys

    :return dict data: Matched data
    """
    data = {}

    def process_value(x: str) -> Union[int, float, str]:
        numerical_value = convert_to_literal(x)
        if numerical_value is None:
            return x.strip()
        else:
            return numerical_value

    try:
        match = re.search(key + '(.+)\n', file_string)
        values = match.group(1).split()
        # Remove escape characters
        parser_key = key.replace('\\', "")
        processed_values = [process_value(raw_value) for raw_value in values]
        data[parser_key] = processed_values[0] if len(processed_values) == 1 else processed_values

    except AttributeError:
        print("Did not find the key:", key)

    if no_colon:
        return {key.rstrip(':'): value for key, value in data.items()}
    else:
        return data


def parse_values_regex(file_string: str, keys: List[str]) -> dict:
    """
    For a list of strings, match the first instance of each string
    (contained in keys) if present in file_string, and return the
    result in a dictionary.

    :param str file_string: Input string
    :param List[str] keys: Keys to match in lines of file_string

    :return dict parsed_data: Matched data
    """
    assert type(keys) == list, "2nd argument of parse_values_regex must be List[str]"
    matches_dict = {}
    for key in keys:
        matches_dict.update(**parse_value_regex(file_string, key))
    return matches_dict
