"""Base class for exciting input classes.
"""
from abc import ABC, abstractmethod
from typing import Union, Set
from xml.etree import ElementTree
from pathlib import Path
import numpy as np

from excitingtools.utils.dict_utils import check_valid_keys


class ExcitingInput(ABC):
    """Base class for exciting inputs."""

    @abstractmethod
    def to_xml(self) -> ElementTree:
        """ Convert class attributes to XML ElementTree."""
        ...


class ExcitingXMLInput(ExcitingInput):
    """Base class for exciting inputs that only consist of many attributes."""

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

    def __init__(self, name: str, valid_attributes: Set[str] = None, **kwargs):
        """Initialise class attributes with kwargs.

        Rather than define all options for a given method, pass as kwargs and directly
        insert as class attributes.

        :param name: Method name.
        """
        self.name = name
        if valid_attributes is not None:
            check_valid_keys(kwargs.keys(), valid_attributes, self.name)
        self.__dict__.update(kwargs)

    def to_xml(self) -> ElementTree:
        """Put class attributes into an XML tree, with the element given by self.name.

        Example ground state XML sub-tree:
           <groundstate vkloff="0.5  0.5  0.5" ngridk="2 2 2" mixer="msec" </groundstate>

        Note, kwargs preserve the order of the arguments, however the order does not appear to be
        preserved when passed to (or perhaps converted to string) with xml.etree.ElementTree.tostring.

        :return ElementTree.Element sub_tree: sub_tree element tree, with class attributes inserted.
        """
        inputs = {}
        attributes = {key: value for key, value in self.__dict__.items() if key != 'name'}
        for key, value in attributes.items():
            inputs[key] = self._attributes_to_input_str[type(value)](value)

        sub_tree = ElementTree.Element(self.name, **inputs)

        # Seems to want this operation on a separate line
        sub_tree.text = ' '

        return sub_tree

    def to_xml_str(self) -> str:
        """ Convert attributes to XML tree string. """
        return ElementTree.tostring(self.to_xml(), encoding='unicode', method='xml')


def query_exciting_version(exciting_root: Union[Path, str]) -> dict:
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
