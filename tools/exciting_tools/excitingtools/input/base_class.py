"""Base class for exciting input classes
"""
from abc import ABC, abstractmethod
from typing import Union
import xml
from xml.etree import ElementTree
from pathlib import Path


class ExcitingInput(ABC):
    """Base class for exciting inputs."""

    @abstractmethod
    def to_xml(self) -> ElementTree:
        ...

    def to_xml_str(self) -> str:
        """ Convert attributes to XML tree string
        """
        return xml.etree.ElementTree.tostring(self.to_xml(), encoding='unicode', method='xml')


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
        raise FileNotFoundError(f'{version_inc} cannot be found. This file generated when the code is built')

    with open(version_inc, 'r') as fid:
        all_lines = fid.readlines()

    git_hash_part1 = all_lines[0].split()[-1][1:-1]
    git_hash_part2 = all_lines[1].split()[-1][1:-1]
    compiler_parts = all_lines[2].split()[2:]
    compiler = " ".join(s for s in compiler_parts).strip()

    version = {'compiler': compiler[1:-1], 'git_hash': git_hash_part1 + git_hash_part2}
    return version
