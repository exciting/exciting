""" Parse exciting species files into dictionary.
"""
from typing import Dict

from excitingtools.parser_utils.parser_decorators import xml_root
from excitingtools.utils.dict_utils import string_value_to_type
from excitingtools.utils.utils import string_to_bool


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
    species = {key: value for key, value in species_tree.attrib.items()}

    for key in ['z', 'mass']:
        species[key] = float(species[key])

    children: Dict[str, list] = {'atomicState': [], 'basis': [], 'muffinTin': []}
    for child in list(species_tree):
        children[child.tag].append(child)

    assert len(children['muffinTin']) == 1, "More than one muffinTin sub-tree in the species file"
    assert len(children['basis']) == 1, "More than one basis sub-tree in the species file"

    muffin_tin_tree = children['muffinTin'][0].attrib
    muffin_tin = {key: float(value) for key, value in muffin_tin_tree.items()}

    atomic_states = []
    for atomic_state_tree in children['atomicState']:
        assert atomic_state_tree.tag == 'atomicState', "Expect tag to be atomicState"
        atomic_states.append(string_value_to_type(atomic_state_tree.attrib))

    basis_tree = children['basis'][0]
    basis: Dict[str, list] = {'default': [], 'custom': [], 'lo': []}

    for func in basis_tree:
        function: dict = func.attrib
        processed_function = string_value_to_type(function)

        if func.tag == 'lo':
            processed_function.update(_parse_lo_from_species(func))

        basis[func.tag].append(processed_function)

    return {
        'species': species,
        'muffin_tin': muffin_tin,
        'atomic_states': atomic_states,
        'basis': basis
        }


def _parse_lo_from_species(lo_function) -> dict:
    """
    Given some lo_function with:
      wf {'matchingOrder': '0', 'trialEnergy': '-2.0', 'searchE': 'true'}
      wf {'matchingOrder': '1', 'trialEnergy': '-2.0', 'searchE': 'true'}

    return
    {'matchingOrder': [0, 1], 'trialEnergy': [-2.0, -2.0], 'searchE': [True, True]}
    """
    # Use lists to GUARANTEE consistent ordering
    matching_order = []
    trial_energy = []
    search = []
    for radial in lo_function:
        matching_order.append(int(radial.attrib.get('matchingOrder')))
        trial_energy.append(float(radial.attrib.get('trialEnergy')))
        search.append(string_to_bool(radial.attrib.get('searchE')))
    return {'matchingOrder': matching_order, 'trialEnergy': trial_energy, 'searchE': search}
