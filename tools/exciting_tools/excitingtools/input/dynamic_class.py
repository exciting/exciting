"""Dynamic Class Construction

class_definitions defined as {'class name':
                                   {'bases': '(ParentClass1, ParentClass2, ...)',
                                    'attributes':
                                             {'attribute1', 'attribute1',
                                              'method1': 'method1'
                                             }
                                    }
                             }

where:
 * 'class name' defines the name of the class, when calling in Python.
 *  'bases' define any parent classes for "class name" to inherit from, given in a tuple.
 *  'attributes' defines any attributes (data or methods) to include in "class name".

ALL data should be defined as a string, as it will get rendered as a Python-valid string,
which subsequently gets interpreted at run-time.
"""
from typing import Dict, Union, List, Optional, Iterator

from excitingtools.exciting_dict_parsers.input_parser import special_tags_to_parse_map
from excitingtools.utils import valid_attributes
from excitingtools.utils.valid_attributes import input_valid_subtrees


def class_constructor_string(class_definitions: Dict[str, Union[dict, str]]) -> str:
    """ Given a dictionary, return a Python-interpretable string of one or more
    class definitions.

        Example function argument:
        class_definitions = {'class name':
                                  {'bases': '(ParentClass1, ParentClass2)',
                                   'attributes':
                                            {'attribute1', 'attribute1',
                                             'method1': 'method1'
                                            }
                                   }
                            }

    :param class_definitions: A dict containing the class name, attributes and parent classes,
    all in string form.
    :return: class_definition_str: A Python-interpretable string of class definitions using
    the type() constructor.
    """
    # first import the base class
    class_definition_str = 'from excitingtools.input.base_class import ExcitingXMLInput \n'
    for class_name, bases_and_attributes in class_definitions.items():
        # Parent class (or classes)
        bases = bases_and_attributes['bases']

        # Attributes and methods (including private attributes)
        attributes = bases_and_attributes['attributes']

        class_definition_str += \
            f"Exciting{class_name}Input = type('Exciting{class_name}Input', {bases}, {attributes}) \n"

    return class_definition_str


def give_class_dictionary(name: str) -> dict:
    """ Gives class dictionary with inheritance and further properties.

    :param name: name of class
    :return: dict with definition
    """
    return {"bases": '(ExcitingXMLInput, )',
            "attributes": {"__doc__": f"Class for exciting {name} input.",
                           "__module__": "excitingtools.input.input_classes",
                           'name': name}}


def class_name_uppercase(name: str) -> str:
    """ Converts the name in a string which better fits as class name, capitalizing the first letter and
    leaving the rest unchanged.
    Exceptions are groundstate -> GroundState, xs -> XS, gw -> GW, bandstructure -> BandStructure, eph -> EPH

    :param name: name of the class/tag in original notation
    :return: name usually capitilized, see exceptions
    """
    exceptions = {"groundstate": "GroundState", "xs": "XS", "gw": "GW", "bandstructure": "BandStructure", "eph": "EPH"}
    if name in exceptions:
        return exceptions[name]
    return name[0].upper() + name[1:]


def get_all_valid_subtrees(valid_xml_tags: Optional[List[str]]) -> Optional[Iterator[str]]:
    """ Get recursively all valid xml (sub-)trees (tags) from exciting.

    :param valid_xml_tags: valid xml tags
    :return: list of all valid xml tags of all trees and subtrees
    """
    if not valid_xml_tags:
        return
    for valid_tag in valid_xml_tags:
        yield valid_tag
        tag_valid_subtrees = valid_attributes.__dict__.get(f"{valid_tag}_valid_subtrees")
        yield from get_all_valid_subtrees(tag_valid_subtrees)


def generate_classes_str() -> str:
    """ Generate string to execute by python. For all standard input classes.

    :return: python string
    """
    non_standard_input_classes = set(special_tags_to_parse_map)
    all_valid_subtrees = set(get_all_valid_subtrees(input_valid_subtrees))
    standard_input_classes = all_valid_subtrees - non_standard_input_classes

    class_definitions = {class_name_uppercase(cls): give_class_dictionary(cls) for cls in standard_input_classes}
    return class_constructor_string(class_definitions)
