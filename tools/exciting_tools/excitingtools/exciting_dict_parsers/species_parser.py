""" Parse exciting species files into dictionary.
"""
from typing import Dict

from excitingtools.parser_utils.parser_decorators import xml_root
from excitingtools.parser_utils.parser_utils import convert_string_dict


@xml_root
def parse_species_xml(root) -> dict:
    """ Parses exciting species files as a dict.

    TODO(Alex) Issue 124. See how easy it is to replace with a generic XML
    parser, with keys defined according to the associated schema.

    Return a dictionary with elements:

      species = {'chemicalSymbol': chemicalSymbol, 'name': name, 'z': z, 'mass': mass}

      muffin_tin = {'rmin': rmin, 'rinf': rinf, 'radius': radius, 'points':  radialmeshPoints}

      atomic_states = [{'n': 1, 'l': 0, 'kappa': 1, 'occ': 2.0, 'core': True},
                      {'n': 2, 'l': 0, 'kappa': 1, 'occ': 2.0, 'core': True}, ...]

      basis['default'] = [{'type': 'lapw', 'trialEnergy': '0.1500', 'searchE': 'true'}]

      basis['custom'] = [{'l': 0, 'type': 'lapw', 'trialEnergy': 1.35670550183736, 'searchE': False},
                         {'l': 1, 'type': 'lapw', 'trialEnergy': -2.69952312512447, 'searchE': False},
                         {'l': 2, 'type': 'lapw', 'trialEnergy': 0.00,  'searchE': False},
                         {'l': 3, 'type': 'lapw', 'trialEnergy': 1.000, 'searchE': False},
                         {'l': 4, 'type': 'lapw', 'trialEnergy': 1.000, 'searchE': False},
                         {'l': 5, 'type': 'lapw', 'trialEnergy': 1.000, 'searchE': False}]

      basis['lo'] = [{'l': 0, 'matchingOrder': [0, 1], 'trialEnergy': [-4.3784, -4.3784], 'searchE': [False, False]},
                     {'l': 0, 'matchingOrder': [0, 1], 'trialEnergy': [1.356, 1.3566], 'searchE': [False, False]},
                    ...]

    :param root: XML file, XML string, or an ET.Element.
    :return : Dictionary of species file data (described above).
    """
    species_tree = root[0]
    species = convert_string_dict(species_tree.attrib)
    
    children: Dict[str, list] = {'atomicState': [], 'basis': [], 'muffinTin': []}
    for child in list(species_tree):
        children[child.tag].append(child)

    assert len(children['muffinTin']) == 1, "More than one muffinTin sub-tree in the species file"
    assert len(children['basis']) == 1, "More than one basis sub-tree in the species file"

    muffin_tin = convert_string_dict(children['muffinTin'][0].attrib)

    atomic_states = []
    for atomic_state_tree in children['atomicState']:
        atomic_states.append(convert_string_dict(atomic_state_tree.attrib))

    basis_tree = children['basis'][0]
    basis: Dict[str, list] = {'default': [], 'custom': [], 'lo': []}

    for func in basis_tree:
        parsed_attributes = convert_string_dict(func.attrib)

        if func.tag == 'lo':
            parsed_attributes["wf"] = [convert_string_dict(wf.attrib) for wf in func]
        
        basis[func.tag].append(parsed_attributes)

    return {
        'species': species,
        'muffin_tin': muffin_tin,
        'atomic_states': atomic_states,
        'basis': basis
        }


