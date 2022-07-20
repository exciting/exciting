""" Parse input XML data directly into corresponding python classes.
"""
from typing import Dict, Union

from excitingtools.input.base_class import ExcitingInput
from excitingtools.input.ground_state import ExcitingGroundStateInput
from excitingtools.input.structure import ExcitingStructure
from excitingtools.input.xs import ExcitingXSInput
from excitingtools.parser_utils.parser_decorators import xml_root
from excitingtools.exciting_dict_parsers.input_parser import root_type, parse_groundstate as parse_groundstate_to_dict, \
    parse_structure as parse_structure_to_dict, parse_xs as parse_xs_to_dict


@xml_root
def parse_groundstate(root: root_type) -> ExcitingGroundStateInput:
    """
    Parse exciting input.xml groundstate element into python ExcitingGroundStateInput.
    :param root: Input for the parser.
    :returns: ExcitingGroundStateInput containing the groundstate input element attributes.
    """
    groundstate: dict = parse_groundstate_to_dict(root)
    return ExcitingGroundStateInput(**groundstate)


@xml_root
def parse_structure(root: root_type) -> ExcitingStructure:
    """
    Parse exciting input.xml structure element into python ExcitingStructure object.
    :param root: Input for the parser.
    :returns: ExcitingStructure containing the structure input element attributes and subelements.
    """
    structure: dict = parse_structure_to_dict(root)
    return ExcitingStructure(
        structure['atoms'],
        structure['lattice'],
        structure['species_path'],
        structure['structure_properties'],
        structure['crystal_properties'],
        structure['species_properties']
    )


@xml_root
def parse_xs(root: root_type) -> Union[ExcitingXSInput, None]:
    """
    Parse exciting input.xml xs element into python ExcitingXSInput object.
    :param root: Input for the parser.
    :returns: ExcitingXSInput containing the xs input element attributes and subelements. Returns None if no xs element
    was found.
    """
    xs: dict = parse_xs_to_dict(root)
    if xs == {}:
        return None
    xs_type = xs.pop('xstype')
    xs_properties = xs.pop('xs_properties')
    return ExcitingXSInput(xs_type, xs_properties, **xs)


exciting_input_type = Union[ExcitingInput, None]


@xml_root
def parse_input_xml(root: root_type) -> Dict[str, exciting_input_type]:
    """
    Parse exciting input.xml into the different python ExcitingInput Objects.
    :param root: Input for the parser.
    :returns: Dictionary which looks like: {'structure': ExcitingStructure,
        'ground_state': ExcitingGroundstateInput, 'xs': ExcitingXSInput}
    If no xs element was found, the value of 'xs' is None.
    """
    structure = parse_structure(root)
    ground_state = parse_groundstate(root)
    xs = parse_xs(root)
    return {'structure': structure, 'groundstate': ground_state, 'xs': xs}
