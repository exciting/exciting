"""Base class for exciting input classes.
"""
import importlib
import warnings
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Union, List, Type, Iterator
from xml.etree import ElementTree

import numpy as np

from excitingtools.exciting_dict_parsers.input_parser import parse_element_xml
from excitingtools.utils import valid_attributes as all_valid_attributes
from excitingtools.utils.dict_utils import check_valid_keys
from excitingtools.utils.jobflow_utils import special_serialization_attrs

path_type = Union[str, Path]


class AbstractExcitingInput(ABC):
    """Base class for exciting inputs."""

    @property
    @abstractmethod
    def name(self) -> str:
        """ Tag of the xml subelement. """
        ...

    @abstractmethod
    def to_xml(self) -> ElementTree:
        """ Convert class attributes to XML ElementTree."""
        ...


class ExcitingXMLInput(AbstractExcitingInput, ABC):
    """Base class for exciting inputs, with exceptions being title, plan and qpointset,
     because they are not passed as a dictionary. """

    # Convert python data to string, formatted specifically for
    _attributes_to_input_str = {int: lambda x: str(x),
                                np.int64: lambda x: str(x),
                                np.float64: lambda x: str(x),
                                float: lambda x: str(x),
                                bool: lambda x: str(x).lower(),
                                str: lambda x: x,
                                list: lambda mylist: " ".join(str(x).lower() for x in mylist).strip(),
                                tuple: lambda mylist: " ".join(str(x).lower() for x in mylist).strip()
                                }

    def __init__(self, **kwargs):
        """Initialise class attributes with kwargs.

        Rather than define all options for a given method, pass as kwargs and directly
        insert as class attributes.

        Valid attributes, subtrees and mandatory attributes are taken automatically from
        the parsed schema, see [valid_attributes.py](excitingtools/utils/valid_attributes.py).
        """
        valid_attributes, valid_subtrees, mandatory_keys = self.get_valid_attributes()

        # check the keys
        missing_mandatory_keys = mandatory_keys - set(kwargs.keys())
        if missing_mandatory_keys:
            raise ValueError(f"Missing mandatory arguments: {missing_mandatory_keys}")
        check_valid_keys(kwargs.keys(), valid_attributes | set(valid_subtrees), self.name)

        # initialise the subtrees
        class_list = self._class_list_excitingtools()
        subtree_class_map = {cls.name: cls for cls in class_list}
        subtrees = set(kwargs.keys()) - valid_attributes
        for subtree in subtrees:
            kwargs[subtree] = self._initialise_subelement_attribute(subtree_class_map[subtree], kwargs[subtree])

        # Set attributes from kwargs
        self.__dict__.update(kwargs)

    def __setattr__(self, name: str, value):
        """ Overload the attribute setting in python with instance.attr = value to check for validity in the schema.

        :param name: name of the attribute
        :param value: new value, can be anything
        """
        valid_attributes, valid_subtrees, _ = self.get_valid_attributes()
        check_valid_keys({name}, valid_attributes | set(valid_subtrees), self.name)
        super().__setattr__(name, value)

    def __delattr__(self, name: str):
        mandatory_keys = list(self.get_valid_attributes())[2]
        if name in mandatory_keys:
            warnings.warn(f"Attempt to delete mandatory attribute '{name}' was prevented.")
        else:
            super().__delattr__(name)

    def get_valid_attributes(self) -> Iterator:
        """ Extract the valid attributes, valid subtrees and mandatory attributes from the parsed schema.

        :return: valid attributes, valid subtrees and mandatory attributes
        """
        yield set(all_valid_attributes.__dict__.get(self.name + "_valid_attributes", set()))
        yield all_valid_attributes.__dict__.get(self.name + "_valid_subtrees", [])
        yield set(all_valid_attributes.__dict__.get(self.name + "_mandatory_attributes", set()))

    @staticmethod
    def _class_list_excitingtools() -> List[Type[AbstractExcitingInput]]:
        """ Find all exciting input classes in own module and excitingtools.
        """
        excitingtools_namespace_content = importlib.import_module("excitingtools").__dict__
        input_class_namespace_content = importlib.import_module("excitingtools.input.input_classes").__dict__
        all_contents = {**excitingtools_namespace_content, **input_class_namespace_content}.values()
        return [cls for cls in all_contents if isinstance(cls, type) and issubclass(cls, AbstractExcitingInput)]

    @staticmethod
    def _initialise_subelement_attribute(XMLClass, element):
        """ Initialize given elements to the ExcitingXSInput constructor. If element is already ExcitingXMLInput class
        object, nothing happens. Else the class constructor of the given XMLClass is called. For a passed
        dictionary the dictionary is passed as kwargs.
        """
        if isinstance(element, XMLClass):
            return element
        elif isinstance(element, dict):
            # assume kwargs
            return XMLClass(**element)
        else:
            # Assume the element type is valid for the class constructor
            return XMLClass(element)

    def to_xml(self) -> ElementTree:
        """Put class attributes into an XML tree, with the element given by self.name.

        Example ground state XML subtree:
           <groundstate vkloff="0.5  0.5  0.5" ngridk="2 2 2" mixer="msec" </groundstate>

        Note, kwargs preserve the order of the arguments, however the order does not appear to be
        preserved when passed to (or perhaps converted to string) with xml.etree.ElementTree.tostring.

        :return ElementTree.Element sub_tree: sub_tree element tree, with class attributes inserted.
        """
        attributes = {key: self._attributes_to_input_str[type(value)](value) for key, value
                      in vars(self).items() if not isinstance(value, AbstractExcitingInput)}
        xml_tree = ElementTree.Element(self.name, **attributes)

        valid_subtrees = all_valid_attributes.__dict__.get(self.name + "_valid_subtrees", [])
        subtrees = {key: self.__dict__[key] for key in set(vars(self).keys()) - set(attributes.keys())}
        ordered_subtrees = [subtrees[x] for x in valid_subtrees if x in subtrees]
        for subtree in ordered_subtrees:
            xml_tree.append(subtree.to_xml())

        # Seems to want this operation on a separate line
        xml_tree.text = ' '

        return xml_tree

    def to_xml_str(self) -> str:
        """ Convert attributes to XML tree string. """
        return ElementTree.tostring(self.to_xml(), encoding='unicode', method='xml')

    def as_dict(self) -> dict:
        """ Convert attributes to dictionary. """
        serialise_attrs = special_serialization_attrs(self)
        return {**serialise_attrs, "xml_string": self.to_xml_str()}

    @classmethod
    def from_xml(cls, xml_string: path_type):
        """ Initialise class instance from XML-formatted string.

        Example Usage
        --------------
        xs_input = ExcitingXSInput.from_xml(xml_string)
        """
        return cls(**parse_element_xml(xml_string, tag=cls.name))

    @classmethod
    def from_dict(cls, d):
        """ Recreates class instance from dictionary. """
        return cls.from_xml(d["xml_string"])


def query_exciting_version(exciting_root: path_type) -> dict:
    """Query the exciting version
    Inspect version.inc, which is constructed at compile-time.

    Assumes version.inc has this structure:
     #define GITHASH "1a2087b0775a87059d53"
     #define GITHASH2 "5d01a5475a10f00d0ad7"
     #define COMPILERVERSION "GNU Fortran (MacPorts gcc9 9.3.0_4) 9.3.0"
     #define VERSIONFROMDATE /21,12,01/

    TODO(Fab) Issue 117. Parse major version.
     Would need to parse src/mod_misc.F90 and regex for "character(40) :: versionname = "
     Refactor whole routine to use regex.

    :param exciting_root: exciting root directory.
    :return version: Build and version details
    """
    if isinstance(exciting_root, str):
        exciting_root = Path(exciting_root)

    version_inc = exciting_root / 'src/version.inc'

    if not version_inc.exists():
        raise FileNotFoundError(f'{version_inc} cannot be found. '
                                f'This file generated when the code is built')

    with open(version_inc, 'r') as fid:
        all_lines = fid.readlines()

    git_hash_part1 = all_lines[0].split()[-1][1:-1]
    git_hash_part2 = all_lines[1].split()[-1][1:-1]
    compiler_parts = all_lines[2].split()[2:]
    compiler = " ".join(s for s in compiler_parts).strip()

    version = {'compiler': compiler[1:-1], 'git_hash': git_hash_part1 + git_hash_part2}
    return version
