"""Utilities to aid in writing and formatting XML
"""
from xml.dom import minidom
from xml.etree import ElementTree
import re


def xml_tree_to_pretty_str(elem: ElementTree.Element) -> str:
    """Convert an XML element to a pretty string.

    :param ElementTree.Element elem: Element/ element tree
    :return str : XML tree string, with pretty formatting.
    """
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="\t")


def line_reformatter(input_str: str, tag: str) -> str:
    """Identify attributes of an XML string element, and reformat them such that they are on new lines.

    NOTE: This could be split into two routines. One that finds the (start, end) of
    each attribute, and one that reformats the string, given these indices.
    (or refactored completely - current solution is not very elegant).

    Example
    ----------
    Takes:
        input_str = <groundstate do="fromscratch" ngridk="6 6 6" nosource="false" rgkmax="8.0"
                     tforce="true" vkloff="0 0 0" xctype="GGA_PBE_SOL"> </groundstate>
        tag = 'groundstate'

    Returns:
        reformatted_str =
        <groundstate
             do="fromscratch"
             ngridk="6 6 6"
             nosource="false"
             rgkmax="8.0"
             tforce="true"
             vkloff="0 0 0"
             xctype="GGA_PBE_SOL">
        </groundstate>

    :param str input_str: Input string, opened and closed with an XML element.
    :param str tag: XML element name.
    :return str reformatted_str: Reformatted form of input_str
    """
    # Get rid of format characters, like \n, \t etc
    input_str = input_str.strip()
    tag = tag.strip()

    # This should not occur if the routine is only used with `prettify_tag_attributes`
    if tag != input_str[0:len(tag)]:
        raise ValueError(f'tag, "{tag}", is inconsistent the element name, "{input_str[0:len(tag)]}"')

    tag_indent = ' ' * 8
    attr_indent = tag_indent + ' ' * 3

    # Isolate attributes according to position of quotation marks in string
    # (cannot use whitespaces to split)
    quote_indices = [x.start() for x in re.finditer('\"', input_str)]
    closing_quote_indices = quote_indices[1::2]
    attribute_start_indices = [len(tag) + 1] + [i + 1 for i in closing_quote_indices[:-1]]

    reformatted_str = tag_indent + tag + '\n'

    for i in range(0, len(attribute_start_indices)):
        i1 = attribute_start_indices[i]
        i2 = closing_quote_indices[i]
        attribute_str = input_str[i1:i2 + 1].strip()
        reformatted_str += attr_indent + attribute_str + '\n'

    reformatted_str += tag_indent + input_str[-(len(tag) + 2):] + '\n'
    return reformatted_str


def prettify_tag_attributes(xml_string: str, tag: str) -> str:
    """Prettify XML string formatting of attributes, for a given XML element.

    The routine finds the line containing the XML element which matches <tag, applies
    a line_reformatter, and replaces the line. If the tag is not matched, the input
    xml_string is returned.

    Example usage:
    ```
        string = <groundstate do="fromscratch" ngridk="6 6 6" nosource="false" rgkmax="8.0" tforce="true" vkloff="0 0 0"
                  xctype="GGA_PBE_SOL"> </groundstate>
        pretty_string = prettify_tag_attributes(string, '<groundstate')
        print(pretty_string)
        > <groundstate
             do="fromscratch"
             ngridk="6 6 6"
             nosource="false"
             rgkmax="8.0"
             tforce="true"
             vkloff="0 0 0"
             xctype="GGA_PBE_SOL">
        </groundstate>
    ```

    :param str xml_string: Already-prettified XML string (assumes tags are on their own lines)
    :param str tag: XML element of the form "<tag"
    :return str reformatted_xml_string: xml_string, with the tag substring reformatted according to the example
     - Line break per attribute.
    """
    xml_list = xml_string.split('\n')

    for i, line in enumerate(xml_list):
        match = re.match(tag, line)
        if match:
            xml_list[i] = line_reformatter(line, tag)
            break

    return "".join(x + '\n' for x in xml_list)
