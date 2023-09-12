""" Tests for utils. """
from excitingtools.utils.utils import can_be_float, convert_to_literal, get_new_line_indices, flatten_list


def test_can_be_float():
    """
    Test can_be_float function on different input types.
    """
    assert can_be_float('1.0'), (
        "Expect string of a literal '1.0' can convert to float")
    assert can_be_float('1'), (
        "Expect string of a literal '1' can convert to float")
    assert can_be_float(True), (
        "Expect True can be converted to a float (would be 1.0)")
    assert can_be_float('a') is False, (
        "Expect string of the letter 'a' cannot be converted to a float")


def test_convert_to_literal():
    """
    Test convert_to_literal turns string reps of numeric data into numeric data.
    """
    assert convert_to_literal('1.1') == 1.1, "string of literal '1.1' converts to float"
    assert convert_to_literal('1.0') == 1.0, "string of literal '1.0' converts to float"
    assert convert_to_literal('1') == 1, "string of literal '1' converts to int"


def test_get_new_line_indices():
    """
    Test getting new line indices function.
    """
    test_string = 'Test Here\n 2nd Line'
    expected_line_indices = 2
    expected_line_index_list = [0, 10]
    assert len(get_new_line_indices(test_string)) == expected_line_indices
    assert get_new_line_indices(test_string)[0] == expected_line_index_list[0]
    assert get_new_line_indices(test_string)[1] == expected_line_index_list[1]


def test_flatten_list():
    input_list = [[1, 2, 3], 4, 5, [6, [7, 8]], {'9': 9, '10': 10}]
    ref_list = [1, 2, 3, 4, 5, 6, 7, 8, {'9': 9, '10': 10}]
    assert list(flatten_list(input_list)) == ref_list
