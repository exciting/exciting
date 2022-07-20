""" Parsers for input.xml.

TODO(Fabian): Issues 117 & 121:
As more sub-elements are implemented in the input files, also add parsers here
"""
import pathlib
import warnings
from typing import Dict, Union
from xml.etree import ElementTree

from excitingtools.parser_utils.parser_decorators import xml_root

# Valid input formats for all parsers
root_type = Union[str, ElementTree.Element, pathlib.Path]


@xml_root
def parse_groundstate(root: root_type) -> dict:
    """
    Parse exciting input.xml groundstate element into python dictionary.
    :param root: Input for the parser.
    :returns: Dictionary containing the groundstate input element attributes.
    """
    ground_state = root.find('groundstate')
    return ground_state.attrib


@xml_root
def parse_structure(root: root_type) -> dict:
    """
    Parse exciting input.xml structure element into python dictionary.
    :param root: Input for the parser.
    :returns: Dictionary containing the structure input element attributes and subelements. Looks like:
        {'atoms': List of atoms with atom positions in fractional coordinates,
         'lattice': List of 3 lattice vectors, 'species_path': species_path as string,
         'structure_properties': dictionary with the structure_properties,
         'crystal_properties': dictionary with the crystal_properties,
         'species_properties': dictionary with the species_properties}
    """
    structure = root.find('structure')
    structure_properties = structure.attrib
    species_path = structure_properties.pop('speciespath')
    crystal = structure.find('crystal')
    crystal_properties = crystal.attrib
    lattice = []
    for base_vect in crystal.findall('basevect'):
        lattice.append([float(x) for x in base_vect.text.split()])

    atoms = []
    species_properties = {}
    for species in structure.findall('species'):
        species_attributes = species.attrib
        species_file = species_attributes.pop('speciesfile')
        species_symbol = species_file[:-4]
        species_properties[species_symbol] = species_attributes
        for atom in species:
            atom_attributes = atom.attrib
            coord = [float(x) for x in atom_attributes.pop('coord').split()]
            atom_dict = {'species': species_symbol, 'position': coord}
            atom_dict.update(atom_attributes)
            atoms.append(atom_dict)

    return {
        'atoms': atoms,
        'lattice': lattice,
        'species_path': species_path,
        'structure_properties': structure_properties,
        'crystal_properties': crystal_properties,
        'species_properties': species_properties
    }


@xml_root
def parse_xs(root: root_type) -> dict:
    """
    Parse exciting input.xml xs element into python dictionary.
    :param root: Input for the parser.
    :returns: Dictionary containing the xs input element attributes and subelements. Could look like:
        {'xstype': xstype as string, 'xs_properties': dictionary with the xs_properties,
         'energywindow': dictionary with the energywindow_properties,
         'screening': dictionary with the screening_properties, 'BSE': dictionary with bse_properties,
         'qpointset': List of qpoints, 'plan': List of tasks}
    """
    xs = root.find('xs')
    if xs is None:
        return {}

    xs_properties = xs.attrib
    xs_type = xs_properties.pop('xstype')
    valid_xml_elements = ['BSE', 'energywindow', 'screening']
    optional_subelements = {}
    for subelement in xs:
        tag = subelement.tag
        if tag == 'qpointset':
            qpointset = []
            for qpoint in subelement:
                qpointset.append([float(x) for x in qpoint.text.split()])
            optional_subelements['qpointset'] = qpointset
        elif tag == 'plan':
            plan = []
            for doonly in subelement:
                plan.append(doonly.attrib.pop('task'))
            optional_subelements['plan'] = plan
        elif tag in valid_xml_elements:
            optional_subelements[tag] = subelement.attrib
        else:
            warnings.warn(f'Subelement {tag} not yet supported. Its ignored...')

    xs_dict = {'xstype': xs_type, 'xs_properties': xs_properties}
    xs_dict.update(optional_subelements)
    return xs_dict


@xml_root
def parse_input_xml(root: root_type) -> Dict[str, dict]:
    """
    Parse exciting input.xml into python dictionaries.
    :param root: Input for the parser.
    :returns: Dictionary which looks like: {'structure': structure_dict,
        'ground_state': groundstate_dict, 'xs': xs_dict}.
    """
    structure = parse_structure(root)
    ground_state = parse_groundstate(root)
    xs = parse_xs(root)
    return {'structure': structure, 'groundstate': ground_state, 'xs': xs}
