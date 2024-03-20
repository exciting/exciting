"""SpeciesFile class.
Provides functionality to read in a species file, add high energy local orbitals (HELOs), and 
write the updated XML structures back to file.
"""
import copy
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Union
from xml.etree import ElementTree
from collections import defaultdict 

from excitingtools.exciting_dict_parsers.species_parser import parse_species_xml
from excitingtools.input.xml_utils import xml_tree_to_pretty_str
from excitingtools.utils.jobflow_utils import special_serialization_attrs


@dataclass
class SpeciesFile:
    """ Holds a species file. """
    species: dict
    muffin_tin: dict
    atomic_states: list
    basis: dict 

    def __post_init__(self):
        """Ensures that the 'basis' dictionary includes 'custom' and 'lo' keys. 
        If these keys are absent, they are initialized as empty lists. 
        """
        self.basis.setdefault('custom', [])
        self.basis.setdefault('lo', [])

    @classmethod
    def from_file(cls, file: Union[str, Path]):
        """ Reads in an exciting species XML file, parses the file as a dict and
        initializes an instance of the class with the parsed data.
        
        :param file: XML species file.
        :return: An instance of the class initialized with the data from the species file.
        """
        return cls(**parse_species_xml(file))
    
    def as_dict(self) -> dict:
        """Returns a dictionary representing the current object for later recreation.
        The serialise attributes are required for recognition by monty and jobflow.
        """
        serialise_attrs = special_serialization_attrs(self)
        return {**serialise_attrs, **self.__dict__}

    @classmethod
    def from_dict(cls, d: dict):
        my_dict = copy.deepcopy(d)
        # Remove key value pairs needed for workflow programs
        # call function on class to get only the keys (values not needed)
        serialise_keys = special_serialization_attrs(cls)
        for key in serialise_keys:
            my_dict.pop(key, None)
        return cls(**my_dict)

    def get_first_helo_n(self, l: int) -> int:
        """Returns the first principle quantum number 'n' for which additional High Energy Local
        Orbitals (HELOs) can be added, for a given angular momentum quantum number 'l'. 
        The 'n' value determination is based on:

        1. For a specific 'l' channel the highest 'n' specified in the atomic states is 'max_n'.
        2. If there exist local orbitals with a 'n' greater than the maximum found in atomic states for 
           the same 'l', 'max_n' is set to this 'n'.
        3. The first added HELO then starts at 'max_n + 1'.
        4. For an 'l' channel not represented in atomic states or as a local orbital, the first added HELO 
            has principle quantum number 'n = l + 1' .

        :param l: angular momentum quantum number 'l' 
        :return: first possible principle quantum number 'n' for a HELO for the given 'l'-channel 
        """
        atomicstate_ns_per_l, lo_ns_per_l = self.get_atomicstate_ns_per_l(), self.get_lo_ns_per_l()
        max_n_for_specified_l = max(atomicstate_ns_per_l.get(l, [l]) + lo_ns_per_l.get(l, [l]))

        # If matchingOrder > 1 for the same l and highest n, skip n+1 HELO.
        ns_per_l_with_mO_over_1 = self.check_matching_orders()
        if ns_per_l_with_mO_over_1.get(l) and max_n_for_specified_l in ns_per_l_with_mO_over_1[l]:
            warnings.warn(f"Warning: HELO skipped for l: {l} and n: {max_n_for_specified_l + 1}")
            return max_n_for_specified_l + 2

        return max_n_for_specified_l + 1
    
    def get_atomicstate_ns_per_l(self) -> Dict[int, list]:
        """Generates a dictionary mapping each 'l' channel present in the atomic states their corresponding 'n' values. 
        Only includes local orbitals with a defined 'n' value. Local orbitals only defined through the attribute 
        'trialEnergy' are ignored.

        :return: Dictionary with 'l' channels as keys and lists of 'n' values from the atomic states.
        """
        atomicstate_ns_per_l = defaultdict(list)
        for state in self.atomic_states:
            atomicstate_ns_per_l[state["l"]].append(state["n"])

        return atomicstate_ns_per_l

    def get_lo_ns_per_l(self) -> Dict[int, list]:
        """Generates a dictionary mapping each 'l' channel present in the local orbitals to their corresponding 'n' 
        values. Only includes local orbitals with a defined 'n' value. Local orbitals only defined through the 
        attribute 'trialEnergy' are ignored.

        :return: Dictionary with 'l' channels as keys and corresponding lists of 'n' values from local orbitals.
        """
        lo_ns_per_l = defaultdict(list)
        for lo in self.basis["lo"]:
            for wf in lo["wf"]:
                if "n" in wf:
                    lo_ns_per_l[lo["l"]].append(wf["n"])

        return lo_ns_per_l
    
    def get_helos_from_species(self) -> Dict[int, set]:
        """Generates a dictionary mapping each 'l' channel to a list of principle quantum numbers 'n' unique to high
        energy/local orbitals (HELOs). This is achieved by comparing 'n' values in local orbitals against those
        in atomic states for each 'l' channel, extracting those exclusive to HELOs.

        :return: Dictionary with 'l' channels as keys and exclusive HELO 'n' values as sets.
        """
        atomicstate_ns_per_l, lo_ns_per_l = self.get_atomicstate_ns_per_l(), self.get_lo_ns_per_l()
        
        helo_ns_per_l = {}
        for l, ns in lo_ns_per_l.items():
            unique_helos = set(ns) - set(atomicstate_ns_per_l.get(l, []))
            helo_ns_per_l[l] = unique_helos

        return helo_ns_per_l

    def check_matching_orders(self) -> Dict[int, list]:
        """Creates a dictionary that lists principle quantum numbers 'n' for each angular momentum quantum number 'l', 
        where the 'n' values are associated with a matchingOrder greater than 1 in the local orbitals.

        :return: Dictionary with 'l' channels as keys and lists of 'n' values that have matchingOrder > 1 as values.
        """
        ns_per_l_with_mO_over_1 = defaultdict(list)
        for lo in self.basis["lo"]:
            for wf in lo["wf"]:
                if wf.get("n") and wf["matchingOrder"] > 1:
                    ns_per_l_with_mO_over_1[lo["l"]].append(wf["n"])

        return ns_per_l_with_mO_over_1

    def add_helos(self, l: int, number: int):
        """Adds a specific 'number' of High Energy Local Orbitals (HELOs) to the basis 
        for a given angular momentum channel 'l'.

        :param l: angular momentum number l
        :param number: the number of HELOs to be added to the l-channel
        """
        first_helo_n_for_l = self.get_first_helo_n(l)

        for nr_lo in range(number):
            new_lo = {
                'l': l,
                'wf': [
                    {'matchingOrder': 0, 'searchE': 'False', 'n': first_helo_n_for_l + nr_lo},
                    {'matchingOrder': 1, 'searchE': 'False', 'n': first_helo_n_for_l + nr_lo}
                ]
            }
            self.basis['lo'].append(new_lo)

    def to_xml(self) -> ElementTree.Element:
        """Converts the class attributes into an XML structure using ElementTree.

        :return: An ElementTree.Element representing the root of the XML structure.
        """
        spdb = ElementTree.Element("spdb")
        sp = ElementTree.SubElement(spdb, "sp", {k: str(v) for k, v in self.species.items()})
        ElementTree.SubElement(sp, "muffinTin", {k: str(v) for k, v in self.muffin_tin.items()})
        for state in self.atomic_states:
            ElementTree.SubElement(
                sp, "atomicState", {k: str(v).lower() for k, v in state.items()}) 

        basis = ElementTree.SubElement(sp, "basis")
        for default in self.basis["default"]:
            ElementTree.SubElement(
                basis, "default", {k: str(v).lower() for k, v in default.items()})

        for custom in self.basis["custom"]:
            ElementTree.SubElement(
                basis, "custom", {k: str(v).lower() for k, v in custom.items()})

        for lo in self.basis["lo"]:
            lo_element = ElementTree.SubElement(basis, "lo", l=str(lo["l"]))
            for wf in lo['wf']:
                ElementTree.SubElement(
                    lo_element, "wf", {k: str(v).lower() for k, v in wf.items()})

        return spdb

    def to_xml_str(self) -> str:
        """Compose XML ElementTrees from exciting input classes to create an input xml string.

        :return: Input XML tree as a string, with pretty formatting.
        """
        return xml_tree_to_pretty_str(self.to_xml())

    def write(self, filename: Union[str, Path]):
        """Writes the xml string to file.

        :param filename: name of the file.
        """
        with open(filename, "w") as fid:
            fid.write(self.to_xml_str())
