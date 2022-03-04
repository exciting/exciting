"""Generate an exciting XML input tree.
"""
from typing import Optional
from collections import OrderedDict
from xml.etree import ElementTree

from excitingtools.input.structure import ExcitingStructure
from excitingtools.input.ground_state import ExcitingGroundStateInput


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
                       title: Optional[str] = '',
                       groundstate: Optional[ExcitingGroundStateInput] = None) -> ElementTree.Element:
    """Compose XML ElementTrees from exciting input classes to create an input XML tree.

    Expected usage: input_xml = exciting_input_xml(structure, title=title, groundstate = groundstate)

    :param ExcitingStructure structure: Structure containing lattice vectors and atomic positions.
    :param Optional[str] title: Optional title for input file
    :param groundstate: exciting ground state input object
    :return ElementTree.Element root: Input XML tree, with sub-elements inserted.
    """
    root = initialise_input_xml(title)

    structure_tree = structure.to_xml()
    root.append(structure_tree)

    exciting_elements = OrderedDict([('groundstate', groundstate)])

    for element in exciting_elements.values():
        if element is not None:
            root.append(element.to_xml())

    return root
