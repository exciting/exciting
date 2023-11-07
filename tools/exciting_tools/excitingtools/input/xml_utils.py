"""Utilities to aid in writing and formatting XML
"""
import re
from xml.dom import minidom
from xml.etree import ElementTree


def xml_tree_to_pretty_str(elem: ElementTree.Element) -> str:
    """Convert an XML element to a pretty string.

    :param ElementTree.Element elem: Element/ element tree
    :return str : XML tree string, with pretty formatting.
    """
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="\t")


def line_reformatter(input_str: str) -> str:
    """Identify attributes of an XML string element, and reformat them such that they are on new lines.

    NOTE: This could be split into two routines. One that finds the (start, end) of
    each attribute, and one that reformats the string, given these indices.
    (or refactored completely - current solution is not very elegant).

    Example
    ----------
    Takes:
        input_str = '\t<groundstate do="fromscratch" ngridk="6 6 6" nosource="false" rgkmax="8.0"
                     tforce="true" vkloff="0 0 0" xctype="GGA_PBE_SOL"> </groundstate>'

    Returns:
        reformatted_str =
        '<groundstate
             do="fromscratch"
             ngridk="6 6 6"
             nosource="false"
             rgkmax="8.0"
             tforce="true"
             vkloff="0 0 0"
             xctype="GGA_PBE_SOL">
        </groundstate>'

    There are 3 possibilities: the xml element has further subelements, than it ends only with a single '>',
    if the xml element has only attributes, than there are two options to close the element: either long with a
    '> </tag>' or short with a '/>'.

    :param str input_str: Input string, opened and closed with an XML element.
    :return str reformatted_str: Reformatted form of input_str
    """
    full_tag = input_str.split(" ")[0]
    tag = full_tag.strip()
    number_of_tag_indents = len(full_tag) - len(tag)

    tag_indent = '\t' * number_of_tag_indents
    attr_indent = tag_indent + ' ' * 3

    # Get rid of format characters, like \n, \t etc
    input_str = input_str.strip()

    # Isolate attributes according to position of quotation marks in string
    # (cannot use whitespaces to split)
    quote_indices = [x.start() for x in re.finditer('\"', input_str)]
    closing_quote_indices = quote_indices[1::2]
    attribute_start_indices = [len(tag) + 1] + [i + 1 for i in
                                                closing_quote_indices[:-1]]

    reformatted_str = tag_indent + tag + '\n'

    for i in range(0, len(attribute_start_indices)):
        i1 = attribute_start_indices[i]
        i2 = closing_quote_indices[i]
        attribute_str = input_str[i1:i2 + 1].strip()
        reformatted_str += attr_indent + attribute_str + '\n'

    # short ending option:
    if input_str[-2:] == '/>':
        return reformatted_str[:-1] + '/>'
    reformatted_str = reformatted_str[:-1] + '>'
    # no ending option:
    if not input_str[-(len(tag) + 2):-len(tag)] == '</':
        return reformatted_str
    # long ending option:
    reformatted_str += "\n" + tag_indent + input_str[-(len(tag) + 2):] + ''
    return reformatted_str


def prettify_tag_attributes(xml_string: str) -> str:
    """Prettify XML string formatting of attributes.

    The routine finds the lines containing an XML element, applies
    a line_reformatter, and replaces the line. Only for long lines with more than 2 attributes.

    Example usage:
    ```
        string = <groundstate do="fromscratch" ngridk="6 6 6" nosource="false" rgkmax="8.0" tforce="true" vkloff="0 0 0"
                  xctype="GGA_PBE_SOL"> </groundstate>
        pretty_string = prettify_tag_attributes(string)
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
    :return str reformatted_xml_string: xml_string, with the tag substrings reformatted according to the example
     - Line break per attribute.
    """
    xml_list = xml_string.split('\n')
    min_number_attributes_for_split = 3

    for i, line in enumerate(xml_list):
        if line.count("=") >= min_number_attributes_for_split:
            xml_list[i] = line_reformatter(line)

    xml_list_without_blank_lines = filter(lambda x: x.strip(), xml_list)
    return "".join(x + '\n' for x in xml_list_without_blank_lines)
