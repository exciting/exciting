""" Parse input XML data directly into corresponding python classes.
"""

from excitingtools.exciting_dict_parsers.input_parser import parse_input_xml as parse_input_xml_to_dict
from excitingtools.input.input_xml import ExcitingInputXML
from excitingtools.parser_utils.parser_decorators import xml_root


@xml_root
def parse_input_xml(root) -> ExcitingInputXML:
    """ Parse exciting input.xml into the corresponding python ExcitingInput Objects.
    :param root: Input for the parser.
    :returns: the exciting input object
    """
    element_dict: dict = parse_input_xml_to_dict(root)
    return ExcitingInputXML(**element_dict)
