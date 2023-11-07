"""Test XML utilities.
"""

from excitingtools.input.xml_utils import line_reformatter


def test_line_reformatter_long_ending():
    """Test reformatting of XML as a string works."""
    input_str = (
        '<groundstate do="fromscratch" ngridk="6 6 6" nosource="false" '
        'rgkmax="8.0" tforce="true" vkloff="0 0 0" xctype="GGA_PBE_SOL"> '
        '</groundstate>'
    )
    pretty_input_str = line_reformatter(input_str)
    # Note, whitespace is important
    reference = """<groundstate
   do="fromscratch"
   ngridk="6 6 6"
   nosource="false"
   rgkmax="8.0"
   tforce="true"
   vkloff="0 0 0"
   xctype="GGA_PBE_SOL">
</groundstate>"""
    assert reference == pretty_input_str


def test_line_reformatter_short_ending():
    """Test reformatting of XML as a string works for xml strings
    with a different (short) ending."""
    input_str = (
        '<groundstate do="fromscratch" ngridk="6 6 6" nosource="false" '
        'rgkmax="8.0" tforce="true" vkloff="0 0 0" xctype="GGA_PBE_SOL"/> '
    )
    pretty_input_str = line_reformatter(input_str)
    # Note, whitespace is important
    reference = """<groundstate
   do="fromscratch"
   ngridk="6 6 6"
   nosource="false"
   rgkmax="8.0"
   tforce="true"
   vkloff="0 0 0"
   xctype="GGA_PBE_SOL"/>"""
    assert reference == pretty_input_str


def test_line_reformatter_no_closing():
    """Test reformatting of XML as a string works for xml strings
    with no closing at the end."""
    input_str = (
        '<groundstate do="fromscratch" ngridk="6 6 6" nosource="false" '
        'rgkmax="8.0" tforce="true" vkloff="0 0 0" xctype="GGA_PBE_SOL"> '
    )
    pretty_input_str = line_reformatter(input_str)
    # Note, whitespace is important
    reference = """<groundstate
   do="fromscratch"
   ngridk="6 6 6"
   nosource="false"
   rgkmax="8.0"
   tforce="true"
   vkloff="0 0 0"
   xctype="GGA_PBE_SOL">"""
    assert reference == pretty_input_str
