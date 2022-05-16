"""Generate an exciting XML input tree.
"""
from typing import Optional
from collections import OrderedDict
from xml.etree import ElementTree

from excitingtools.input.structure import ExcitingStructure
from excitingtools.input.ground_state import ExcitingGroundStateInput
from excitingtools.input.xs import ExcitingXSInput
from excitingtools.input.xml_utils import xml_tree_to_pretty_str, prettify_tag_attributes


def initialise_input_xml(title: str) -> ElementTree.Element:
    """Initialise input.xml element tree for exciting.

    Information on the schema reference ignored, but could be reintroduced for validation purposes.
    root.set(
       '{http://www.w3.org/2001/XMLSchema-instance}noNamespaceSchemaLocation',
       'http://xml.exciting-code.org/excitinginput.xsd')

    :param str title: Title for calculation.
    :return ElementTree.Elementroot: Element tree root.
    """
    root = ElementTree.Element('input')
    ElementTree.SubElement(root, 'title').text = title
    return root


def exciting_input_xml(structure: ExcitingStructure,
                       groundstate: ExcitingGroundStateInput,
                       title: Optional[str] = '',
                       xs: Optional[ExcitingXSInput] = None) -> ElementTree.Element:
    """Compose XML ElementTrees from exciting input classes to create an input XML tree.

    Expected usage: input_xml = exciting_input_xml(structure, groundstate, title=title)

    :param ExcitingStructure structure: Structure containing lattice vectors and atomic positions.
    :param groundstate: exciting ground state input object.
    :param Optional[str] title: Optional title for input file.
    :param xs: exciting xs input object.
    :return ElementTree.Element root: Input XML tree, with sub-elements inserted.
    """
    root = initialise_input_xml(title)

    structure_tree = structure.to_xml()
    root.append(structure_tree)

    exciting_elements = OrderedDict([('groundstate', groundstate), ('xs', xs)])

    for element in exciting_elements.values():
        if element is not None:
            root.append(element.to_xml())

    return root


def exciting_input_xml_str(structure: ExcitingStructure,
                           groundstate: ExcitingGroundStateInput,
                           **kwargs) -> str:
    """ "Compose XML ElementTrees from exciting input classes to create an input xml string.

    :param ExcitingStructure structure: Structure containing lattice vectors and atomic positions.
    :param groundstate: exciting ground state input object.
    :return input_xml_str: Input XML tree as a string, with pretty formatting.
    """
    xml_tree = exciting_input_xml(structure, groundstate, **kwargs)
    tags_to_prettify = ["\t<structure", "\t\t<crystal", "\t\t<species", "\t\t\t<atom", "\t<groundstate", "\t<xs",
                        "\t\t<screening", "\t\t<BSE", "\t\t<energywindow"]
    input_xml_str = prettify_tag_attributes(xml_tree_to_pretty_str(xml_tree), tags_to_prettify)
    return input_xml_str
