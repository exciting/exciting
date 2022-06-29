"""Test XML utilities.
"""
import pytest

from excitingtools.input.xml_utils import line_reformatter, \
    prettify_tag_attributes


def test_line_reformatter_long_ending():
    """Test reformatting of XML as a string works."""
    input_str = (
        '<groundstate do="fromscratch" ngridk="6 6 6" nosource="false" '
        'rgkmax="8.0" tforce="true" vkloff="0 0 0" xctype="GGA_PBE_SOL"> '
        '</groundstate>'
    )
    pretty_input_str = line_reformatter(input_str, '<groundstate')
    # Note, whitespace is important
    reference = """<groundstate
   do="fromscratch"
   ngridk="6 6 6"
   nosource="false"
   rgkmax="8.0"
   tforce="true"
   vkloff="0 0 0"
   xctype="GGA_PBE_SOL">
</groundstate>
"""
    assert reference == pretty_input_str


def test_line_reformatter_short_ending():
    """Test reformatting of XML as a string works for xml strings
    with a different (short) ending."""
    input_str = (
        '<groundstate do="fromscratch" ngridk="6 6 6" nosource="false" '
        'rgkmax="8.0" tforce="true" vkloff="0 0 0" xctype="GGA_PBE_SOL"/> '
    )
    pretty_input_str = line_reformatter(input_str, '<groundstate')
    # Note, whitespace is important
    reference = """<groundstate
   do="fromscratch"
   ngridk="6 6 6"
   nosource="false"
   rgkmax="8.0"
   tforce="true"
   vkloff="0 0 0"
   xctype="GGA_PBE_SOL"/>
"""
    assert reference == pretty_input_str


def test_line_reformatter_no_closing():
    """Test reformatting of XML as a string works for xml strings
    with no closing at the end."""
    input_str = (
        '<groundstate do="fromscratch" ngridk="6 6 6" nosource="false" '
        'rgkmax="8.0" tforce="true" vkloff="0 0 0" xctype="GGA_PBE_SOL"> '
    )
    pretty_input_str = line_reformatter(input_str, '<groundstate')
    # Note, whitespace is important
    reference = """<groundstate
   do="fromscratch"
   ngridk="6 6 6"
   nosource="false"
   rgkmax="8.0"
   tforce="true"
   vkloff="0 0 0"
   xctype="GGA_PBE_SOL">
"""
    assert reference == pretty_input_str


def test_line_reformatter_unmatched_tag():
    """
    Test error occurs when tag is inconsistent with the element name.

    Recall the tag is passed to `line_reformatter`. This should not occur if
    `line_reformatter` is only called from `prettify_tag_attributes`.
    """
    input_str = (
        '<groundstate do="fromscratch" ngridk="6 6 6" nosource="false" '
        'rgkmax="8.0" tforce="true" vkloff="0 0 0" xctype="GGA_PBE_SOL"> '
        '</groundstate>'
    )

    with pytest.raises(ValueError) as error:
        line_reformatter(input_str, '<error')

    assert error.value.args[0] == 'tag, "<error", is inconsistent the element name, "<groun"'


def test_prettify_tag_attributes():
    """Test prettifying the tag attributes."""
    input_str = (
        '<groundstate do="fromscratch" ngridk="6 6 6" nosource="false" '
        'rgkmax="8.0" tforce="true" vkloff="0 0 0" xctype="GGA_PBE_SOL"> '
        '</groundstate>'
    )

    output = prettify_tag_attributes(input_str, '<unmatched')
    assert output.rstrip() == input_str, (
        "If the tag is not matched in the XML string, the input string should "
        "be returned - (newline character stripped from end)")
