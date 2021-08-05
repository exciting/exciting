from excitingtools.utils import can_be_float, convert_to_literal


def test_can_be_float():
    assert can_be_float('1.0'), "Expect string of a literal '1.0' can convert to float"
    assert can_be_float('1'), "Expect string of a literal '1' can convert to float"
    assert can_be_float(True), "Expect True can be converted to a float (would be 1.0)"
    assert can_be_float('a') == False, "Expect string of the letter 'a' cannot be converted to a float"


def test_convert_to_literal():
    assert convert_to_literal('1.1') == 1.1, "string of literal '1.1' converts to float"
    assert convert_to_literal('1.0') == 1.0, "string of literal '1.0' converts to float"
    assert convert_to_literal('1') == 1, "string of literal '1' converts to int"
