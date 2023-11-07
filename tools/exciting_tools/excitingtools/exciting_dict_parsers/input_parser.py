""" Parsers for input.xml.

TODO(Fabian): Issues 117 & 121:
As more sub-elements are implemented in the input files, also add parsers here
"""
import copy
from typing import Tuple
from xml.etree import ElementTree

from excitingtools.parser_utils.parser_decorators import xml_root
from excitingtools.parser_utils.parser_utils import find_element, convert_string_dict
from excitingtools.utils import valid_attributes as all_valid_attributes
from excitingtools.utils.valid_attributes import input_valid_attributes


@xml_root
def parse_input_xml(root):
    """ Parse an input.xml file into dictionary. """
    assert root.tag == "input"
    return parse_element_xml(root)


def get_root_from_tag(root: ElementTree.Element, tag: str = None) -> Tuple[ElementTree.Element, str]:
    """ Get the root from a tag.

    :param tag: tag of interest
    :param root: xml root containing the tag (or having the specified tag as tag)
    :returns: the tag and the found root, if tag was None returns the tag of the given root
    """
    if tag is None:
        return root, root.tag

    root = find_element(root, tag)
    if root is None:
        raise ValueError(f"Your specified input has no tag {tag}.")

    return root, root.tag


@xml_root
def parse_element_xml(root, tag: str = None) -> dict:
    """ Parse a xml element into dictionary. Can be input.xml root or a subelement of it.
    Put the attributes simply in dict and add recursively the subtrees and nested dicts.

    :param tag: the tag to parse
    :param root: the xml root containing the tag
    :returns: the parsed dictionary, data converted to actual data types
    """
    root, tag = get_root_from_tag(root, tag)

    if tag in special_tags_to_parse_map.keys():
        return special_tags_to_parse_map[tag](root)

    element_dict = convert_string_dict(copy.deepcopy(root.attrib))

    multiple_children = set(getattr(all_valid_attributes, f"{tag}_multiple_children", set()))
    subelements = list(root)
    multiple_tags = set([subelement.tag for subelement in subelements]) & multiple_children
    element_dict.update({tag: [] for tag in multiple_tags})

    # for loop over all subelements to retain the order
    for subelement in subelements:
        if subelement.tag in multiple_children:
            element_dict[subelement.tag].append(parse_element_xml(subelement))
        else:
            element_dict[subelement.tag] = parse_element_xml(subelement)

    return element_dict


def _parse_input_tag(root) -> dict:
    """ Parse special input tag. Necessary because exciting/xml allows arbitrary attributes at this level.
    Only parses the explicitly named attributes in the schema.

    :param root: the xml root containing the input tag
    :returns: the parsed dictionary, data converted to actual data types
    """
    valid_attribs = {key: value for key, value in root.attrib.items() if key in input_valid_attributes}
    element_dict = convert_string_dict(valid_attribs)

    subelements = list(root)
    for subelement in subelements:
        element_dict[subelement.tag] = parse_element_xml(subelement)

    return element_dict


@xml_root
def parse_structure(root) -> dict:
    """ Parse exciting input.xml structure element into python dictionary.

    :param root: Input for the parser.
    :returns: Dictionary containing the structure input element attributes and subelements. Looks like:
        {'atoms': List of atoms with atom positions in fractional coordinates,
         'lattice': List of 3 lattice vectors, 'species_path': species_path as string,
         'crystal_properties': dictionary with the crystal_properties,
         'species_properties': dictionary with the species_properties,
         all additional keys are structure attributes}
    """
    structure = find_element(root, 'structure')
    structure_properties = convert_string_dict(copy.deepcopy(structure.attrib))
    species_path = structure_properties.pop('speciespath')

    crystal = structure.find('crystal')
    crystal_properties = convert_string_dict(copy.deepcopy(crystal.attrib))
    lattice = []
    for base_vect in crystal:
        lattice.append([float(x) for x in base_vect.text.split()])

    atoms = []
    species_properties = {}
    for species in structure.findall('species'):
        species_attributes = convert_string_dict(copy.deepcopy(species.attrib))
        species_file = species_attributes.pop('speciesfile')
        species_symbol = species_file[:-4]

        species_subelements = list(species)
        atom_xml_trees = [x for x in species_subelements if x.tag == "atom"]
        for atom in atom_xml_trees:
            atom_attributes = convert_string_dict(copy.deepcopy(atom.attrib))
            atom_dict = {'species': species_symbol, 'position': atom_attributes.pop('coord')}
            atom_dict.update(atom_attributes)
            atoms.append(atom_dict)

        other_xml_trees = set(species_subelements) - set(atom_xml_trees)
        for tree in other_xml_trees:
            species_attributes[tree.tag] = parse_element_xml(tree)
        species_properties[species_symbol] = species_attributes

    return {
        'atoms': atoms,
        'lattice': lattice,
        'species_path': species_path,
        'crystal_properties': crystal_properties,
        'species_properties': species_properties,
        **structure_properties
    }


# special tag to parse function map or lambda if one-liner
# necessary for tags which doesn't contain simply xml attributes and subtrees
special_tags_to_parse_map = {"input": _parse_input_tag,
                             "title": lambda root: root.text,
                             "structure": parse_structure,
                             "qpointset": lambda root: [[float(x) for x in qpoint.text.split()] for qpoint in root],
                             "plan": lambda root: [doonly.attrib['task'] for doonly in root]}
