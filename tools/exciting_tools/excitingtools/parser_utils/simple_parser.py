""" Simple text parsers.
"""
from typing import Union


def match_current_return_line_n(file_string: str, match: str, n_line=1) -> Union[str, None]:
    """
    Match a string on the ith line and return a substring from the i+n_line line.

    :param str file_string: Input string
    :param str match: String to match
    :param int n_line: The index of the line to return, following the matched line

    :return Union[str, None] matched line string, or None
    """
    file = file_string.split('\n')
    for i, line in enumerate(file):
        if match in line:
            return file[i + n_line]
    return None


def match_current_extract_from_line_n(input_string: str,
                                      keys_extractions: dict,
                                      n_line=1) -> dict:
    """
    Given an input_string, match a substring (defined by the key of keys_extractions),
    return the a substring from the i+n_line line below the match, and extract a value
    from that line using the value of the key, defined by keys_extractions[key]

    :param str file_string: Input string
    :param dict keys_extractions: keys = strings to match
                                  values = methods of data extraction from the strings
    :param int n_line: The index of the line to return, following the matched line

    :return dict data: Dictionary of extracted data from the nth line below a matched string.
    """
    data = {}
    for key, extract_data in keys_extractions.items():
        match = match_current_return_line_n(input_string, key, n_line=n_line)
        if match is not None:
            data[key] = extract_data(match)
    return data
