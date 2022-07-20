"""
Test parsers in simple_parser
"""

from excitingtools.parser_utils.simple_parser import match_current_return_line_n, \
    match_current_extract_from_line_n


def test_match_current_return_next_n():

    input_string = """
       Tolerance factor to reduce the MB size based on
       the eigenvectors of the bare Coulomb potential:   0.100000000000000     
     
    --------------------------------------------------------------------------------
     
     Screened Coulomb potential:
       Full-frequency Random-Phase Approximation
       Some additional data : 5
     """

    string_match = match_current_return_line_n(
        input_string, 'Screened Coulomb potential:', n_line=1
        )
    assert string_match.strip() == "Full-frequency Random-Phase Approximation", \
        "Expect to return the string one line below the match"

    string_match = match_current_return_line_n(input_string, 'Screened Coulomb potential:')
    assert string_match.strip() == "Full-frequency Random-Phase Approximation", \
        "Expect default behaviour to return the string one line below the match"

    string_match = match_current_return_line_n(
        input_string, 'Screened Coulomb potential:', n_line=2
        )
    assert string_match.strip() == "Some additional data : 5", \
        "Expect to return the string two lines below the match"


def test_match_current_extract_from_line_n():

    input_string = """
       Tolerance factor to reduce the MB size based on
       the eigenvectors of the bare Coulomb potential:   0.100000000000000     

    --------------------------------------------------------------------------------

     Screened Coulomb potential:
       Full-frequency Random-Phase Approximation
       Some additional data : 5
     """

    # Extract 5 from "Some additional data : 5"
    key_extraction = {'Screened Coulomb potential:': lambda x: int(x.split(':')[-1])}

    data = match_current_extract_from_line_n(input_string, key_extraction, n_line=2)

    assert data == {'Screened Coulomb potential:': 5}, \
        "Expect to return the integer from 2 lines below the match key"
