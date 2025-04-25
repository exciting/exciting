"""SpeciesFile class.
Provides functionality to read in a species file, add high energy local orbitals (HELOs), and
write the updated XML structures back to file.
"""

import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Set, Tuple, Union
from xml.etree import ElementTree

from excitingtools.base import ECTObject
from excitingtools.exciting_dict_parsers.species_parser import parse_species_xml
from excitingtools.input.xml_utils import xml_tree_to_pretty_str


@dataclass(frozen=True)
class LocalOrbital(ECTObject):
    """A local orbital. Defined by a state and a matching order."""

    species: str
    l_value: int
    n_values: Tuple[int, int]
    matching_orders: Tuple[int, int]


LocalOrbitals = List[LocalOrbital]

QuantumNumberMapping = Dict[int, set]


@dataclass
class SpeciesFile(ECTObject):
    """Holds a species file."""

    species: dict
    muffin_tin: dict
    atomic_states: list
    basis: dict

    def __post_init__(self):
        """Ensures that the 'basis' dictionary includes 'custom' and 'lo' keys.
        If these keys are absent, they are initialized as empty lists.
        """
        self.basis.setdefault("custom", [])
        self.basis.setdefault("lo", [])

    @classmethod
    def from_file(cls, file: Union[str, Path]):
        """Reads in an exciting species XML file, parses the file as a dict and
        initializes an instance of the class with the parsed data.

        :param file: XML species file.
        :return: An instance of the class initialized with the data from the species file.
        """
        return cls(**parse_species_xml(file))

    def get_first_helo_n(self, l: int, skip_lin_dep: bool = True) -> int:
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
        :param skip_lin_dep: If set to True and matchingOrder > 1 for HELO l and n, skip n+1 HELO
        :return: first possible principle quantum number 'n' for a HELO for the given 'l'-channel
        """
        atomicstate_ns_per_l, lo_ns_per_l = self.get_atomicstates_ns_per_l(), self.get_lo_ns_per_l()
        max_n_for_specified_l = max(atomicstate_ns_per_l.get(l, {l}) | lo_ns_per_l.get(l, {l}))

        # If matchingOrder > 1 for the same l and highest n, skip n+1 HELO.
        ns_per_l_with_mO_over_1 = self.check_matching_orders()
        if skip_lin_dep and ns_per_l_with_mO_over_1.get(l) and max_n_for_specified_l in ns_per_l_with_mO_over_1[l]:
            logging.info(f"HELO skipped for l: {l} and n: {max_n_for_specified_l + 1}")
            return max_n_for_specified_l + 2

        return max_n_for_specified_l + 1

    def get_atomicstates_ns_per_l(self, cond: Callable[[dict], bool] = lambda _: True) -> QuantumNumberMapping:
        """Generates a dictionary mapping each 'l' channel present in the atomic states their corresponding 'n' values.

        :param cond: condition for choosing specific atomic states
        :return: Dictionary with 'l' channels as keys and sets of 'n' values from the atomic states.
        """
        atomicstate_ns_per_l = defaultdict(set)
        cond_atomic_states = filter(cond, self.atomic_states)
        for state in cond_atomic_states:
            atomicstate_ns_per_l[state["l"]].add(state["n"])

        return atomicstate_ns_per_l

    def get_valence_and_semicore_atomicstate_ns_per_l(self) -> Dict[int, Dict[str, Set[int]]]:
        """
        Classifies and returns all non-core states into valence and semicore states for each 'l' channel and
        corresponding 'n' value.

        This function classifies all states that do not belong to the core based on the maximum 'n' value for each 'l':
            - States with the maximum 'n' value for their 'l' channel are classified as valence states.
            - All other non-core states are classified as semicore states.

        :return: A dictionary mapping each 'l' channel to another dictionary containing two keys: 'semicore' and
                'valence'. Each of these keys maps to a set of 'n' values corresponding to the classified atomic states.

        Example:
            {
                0: {'semicore': {1, 2}, 'valence': {3}},
                1: {'semicore': {2}, 'valence': {3, 4}},
                ...
            }
        """
        noncore_atomicstate_ns_per_l = self.get_atomicstates_ns_per_l(lambda x: not x["core"])
        valence_and_semicore_states = {}

        for l, ns in noncore_atomicstate_ns_per_l.items():
            max_n = max(ns)
            ns.remove(max_n)
            valence_and_semicore_states[l] = {"semicore": ns, "valence": {max_n}}

        return valence_and_semicore_states

    def get_lo_ns_per_l(self) -> QuantumNumberMapping:
        """Generates a dictionary mapping each 'l' channel present in the local orbitals to their corresponding 'n'
        values. Only includes local orbitals with a defined 'n' value. Local orbitals only defined through the
        attribute 'trialEnergy' are ignored.

        :return: Dictionary with 'l' channels as keys and corresponding sets of 'n' values from local orbitals.
        """
        lo_ns_per_l = defaultdict(set)
        for lo in self.basis["lo"]:
            for wf in lo["wf"]:
                if "n" in wf:
                    lo_ns_per_l[lo["l"]].add(wf["n"])

        return lo_ns_per_l

    def get_los_from_species_file(self) -> LocalOrbitals:
        """
        Generates a list of local orbitals in the species file. Each local orbital is represented by a set of the form
        (l, [n1, n2], [mO1, mO2]). Local orbitals added to the species file in a custom apw+lo element are not
        cosidered.

        :return: List of local orbitals
        """
        los = []
        for lo in self.basis["lo"]:
            ns, mOs = [], []
            for wf in lo["wf"]:
                ns.append(wf["n"])
                mOs.append(wf["matchingOrder"])
            los.append(LocalOrbital(self.species["chemicalSymbol"], lo["l"], tuple(ns), tuple(mOs)))

        return los

    def get_custom_los_from_species_file(self) -> LocalOrbitals:
        """
        Generates a list of local orbitals in the species file as a custom apw+lo element.
        Each local orbital is represented by a set of the form (l, [n1, n2], [mO1, mO2]).

        :return: List of local orbitals
        """

        return [
            LocalOrbital(self.species["chemicalSymbol"], custom["l"], (custom["n"], custom["n"]), (0, 1))
            for custom in self.basis["custom"]
            if custom["type"] == "apw+lo"
        ]

    def get_helos_from_species(self) -> QuantumNumberMapping:
        """Generates a dictionary mapping each 'l' channel to a list of principle quantum numbers 'n' unique to high
        energy/local orbitals (HELOs). This is achieved by comparing 'n' values in local orbitals against those
        in atomic states for each 'l' channel, extracting those exclusive to HELOs.

        :return: Dictionary with 'l' channels as keys and exclusive HELO 'n' values as sets.
        """
        atomicstate_ns_per_l, lo_ns_per_l = self.get_atomicstates_ns_per_l(), self.get_lo_ns_per_l()

        helo_ns_per_l = {}
        for l, ns in lo_ns_per_l.items():
            unique_helos = ns - atomicstate_ns_per_l.get(l, set())
            helo_ns_per_l[l] = unique_helos

        return helo_ns_per_l

    def check_matching_orders(self) -> QuantumNumberMapping:
        """Creates a dictionary that lists principle quantum numbers 'n' for each angular momentum quantum number 'l',
        where the 'n' values are associated with a matchingOrder greater than 1 in the local orbitals.

        :return: Dictionary with 'l' channels as keys and lists of 'n' values that have matchingOrder > 1 as values.
        """
        ns_per_l_with_mO_over_1 = defaultdict(set)

        for lo in self.basis["lo"]:
            for wf in lo["wf"]:
                if wf.get("n") and wf["matchingOrder"] > 1:
                    ns_per_l_with_mO_over_1[lo["l"]].add(wf["n"])

        return ns_per_l_with_mO_over_1

    def add_number_los_for_all_valence_semicore_states(self, number_lo: int):
        """Adds a certain number `number_lo` of local orbitals with increasing matching order for every valence and
        semicore state to the basis `self.basis[lo]`. Here the local orbital always consists of 2 wave function
        elements.

        :param number_lo: number of local orbitals added for every valence and semicore state
        """
        valence_semicore_states = self.get_atomicstates_ns_per_l(lambda x: not x["core"])

        for l, ns in valence_semicore_states.items():
            for n in ns:
                for _ in range(number_lo):
                    self.add_lo_higher_matching_order(l, n)

    def add_basic_lo_all_semicore_states(self):
        """Adds one local orbital for every semicore state present in the atomic states.

        Moving a state from core to valence requires the addition of at least one local orbital to accurately
        describe the semicore state.
        The local orbital consists of two radial functions at matchingOrder 0:
            - one at the same principal quantum number 'n' of the (L)APW  for the given angular momentum ('l') channel.
            - the other at the principal quantum number 'n' of the semi core state.
        """
        valence_and_semicore_states = self.get_valence_and_semicore_atomicstate_ns_per_l()
        semicore_states = {l: data["semicore"] for l, data in valence_and_semicore_states.items()}

        for l, ns in semicore_states.items():
            for n in ns:
                self.add_lo(l, (n + 1, n), (0, 0))

    def add_custom_for_all_valence_states(self, custom_type: str):
        """Adds a custom element to the basis for every valence state present in the atomic states.

        :param custom_type: Type of the custom basis function. It can be either `apw`, `lapw`, or `apw+lo`.
        """
        valence_and_semicore_states = self.get_valence_and_semicore_atomicstate_ns_per_l()
        valence_states = {l: data["valence"] for l, data in valence_and_semicore_states.items()}

        for l, ns in valence_states.items():
            for n in ns:
                self.basis["custom"].append({"l": l, "type": custom_type, "n": n, "searchE": False})

    def add_helos(self, l: int, number: int, skip_lin_dep: bool = True):
        """Adds a specific 'number' of High Energy Local Orbitals (HELOs) to the basis
        for a given angular momentum channel 'l'.

        :param l: angular momentum number l
        :param number: the number of HELOs to be added to the l-channel
        :param skip_lin_dep: If set to True and matchingOrder > 1 for HELO l and n, skip n+1 HELO
        """
        first_helo_n_for_l = self.get_first_helo_n(l, skip_lin_dep)

        for nr_lo in range(number):
            self.add_lo(l, (first_helo_n_for_l + nr_lo,) * 2, [0, 1])

    def find_highest_matching_order_for_state(self, l: int, n: int) -> int:
        """Returns the highest matching order for a specific state, defined by principal quantum numner 'n' and
        anguar momentum 'l', for which a local orbital exists in the species file.

        :param l: angular momentum number l
        :param n: principal quantum number n
        :return: highest matchingOrder mO
        """

        l_los = filter(lambda x: x["l"] == l, self.basis["lo"])
        mOs = [wf["matchingOrder"] for lo in l_los for wf in lo["wf"] if wf.get("n") == n]

        return max(mOs, default=0)

    def add_lo_higher_matching_order(self, l: int, n: int, *, raise_exception: bool = True) -> None:
        """Adds a Local Orbital with the next highest matching order for a state defined by angular momentum 'l'
        and principal quantum number 'n'.

        :param l: angular momentum number
        :param n: principal quantum number
        :param raise_exception: if true, raises an expection if maximum matching order reached
        """
        mO = self.find_highest_matching_order_for_state(l, n)
        self.add_lo(l, (n, n), (mO, mO + 1), raise_exception=raise_exception)

    def add_lo(
        self, l: int, ns: Tuple[int, int], matching_orders: Tuple[int, int], *, raise_exception: bool = True
    ) -> None:
        """Adds a single local orbital to the basis for a given angular momentum channel 'l',
        with tuples of principal quantum number 'n' and corresponding 'matching_orders'.

        :param l: angular momentum number l
        :param ns: tuple of principal quantum number n
        :param matching_orders: tuple of matching orders
        :param raise_exception: if true, raises an expection if maximum matching order reached
        """
        assert len(ns) == len(matching_orders), (
            "Number of principal quantum numbers n must equal the number of given matching orders."
        )

        if matching_orders[0] >= 3:
            if raise_exception:
                raise ValueError("Maximum matchingOrder reached; cannot add new local orbital.")
            return

        wf = [{"matchingOrder": mO, "searchE": False, "n": n} for mO, n in zip(matching_orders, ns)]
        self.basis["lo"].append({"l": l, "wf": wf})

    def remove_lo(
        self, l: int, ns: Tuple[int, int], matching_orders: Tuple[int, int], *, raise_exception: bool = True
    ) -> None:
        """Removes a single local orbital from the basis for a given angular momentum channel 'l',
        based on a list of tuples, each containing a principal quantum number 'n' and its corresponding matching order.
        If a local orbital with the specified 'n' and matching orders is present in the basis, it is removed.

        :param l: Angular momentum number l.
        :param ns: tuple of principal quantum number n
        :param matching_orders: tuple of matching orders
        :param raise_exception: if true, raises an expection if it is not possible to remove the local orbital
        """
        for lo in self.basis["lo"]:
            ns_mOs_lo = sorted([(wf["n"], wf["matchingOrder"]) for wf in lo["wf"] if wf.get("n")])
            if lo["l"] == l and sorted(zip(ns, matching_orders)) == ns_mOs_lo:
                self.basis["lo"].remove(lo)
                return

        if raise_exception:
            raise ValueError("Could not remove local orbital.")

    def add_default(self, trial_energy: float, default_type: str):
        """Adds the default element with a given trial energy and type to the basis.

        :param trial_energy: trial energy used in the default element
        :param default_type: type of the default basis functions. Can be either `apw`, `lapw`, or `apw+lo`.
        """
        self.basis["default"].append({"type": default_type, "trialEnergy": trial_energy, "searchE": False})

    def add_custom_for_high_l(self, highest_valence_l: int, n_high_l: int) -> None:
        """Add custom elements of type "LAPW" for a number of unoccupied l-channels to species file.

        :param highest_valence_l: the highest valence l value
        :param n_high_l: number of unoccupied l-channels for which custom LAPWs are added
        """
        for i in range(n_high_l):
            new_l = highest_valence_l + 1 + i
            self.basis["custom"].append({"l": new_l, "type": "lapw", "n": new_l + 1, "searchE": False})

    def add_lo_for_high_l(
        self, highest_valence_l: int, n_high_l: int, n_high_n: int, highest_mO: int, skip_lin_dep: bool = True
    ) -> None:
        """Add local orbitals for a number of unoccupied l-channels and n-channels. They are added up to a highest
        matching order.

        :param highest_valence_l: the highest valence l value
        :param n_high_l: number of unoccupied l-channels for which local orbitals are added
        :param n_high_n: number of n values per unoccupied l-channel for which local orbitals are added
        :param highest_mO: highest matching order up to which local orbitals are added.
        :param skip_lin_dep: If set to True and matchingOrder > 1 for HELO l and n, skip n+1 HELO
        """
        for li in range(n_high_l):
            l_val = highest_valence_l + li + 1
            first_helo_n_for_l = self.get_first_helo_n(l_val, skip_lin_dep)
            for nr_lo in range(n_high_n):
                if (first_helo_n_for_l + nr_lo) != l_val + 1:
                    self.add_lo(l_val, (first_helo_n_for_l + nr_lo - 1, first_helo_n_for_l + nr_lo), (0, 0))
                for mO in range(highest_mO):
                    self.add_lo(l_val, (first_helo_n_for_l + nr_lo,) * 2, (mO, mO + 1))

    def add_high_n_lo_for_all_valence_states(
        self, highest_valence_l: int, n_high_n: int, highest_mO: int, skip_lin_dep: bool = True
    ) -> None:
        """Add high-n local orbitals for all occupied l-channels. They are added up to a highest
        matching order.

        :param highest_valence_l: the highest valence l_val value
        :param n_high_n: number of n values (above valence n) for which local orbitals are added.
        :param highest_mO: highest matching order up to which local orbitals are added.
        :param skip_lin_dep: If set to True and matchingOrder > 1 for HELO l and n, skip n+1 HELO
        """
        for l_val in range(highest_valence_l + 1):
            first_helo_n_for_l = self.get_first_helo_n(l_val, skip_lin_dep)
            for nr_lo in range(n_high_n):
                self.add_lo(l_val, (first_helo_n_for_l + nr_lo - 1, first_helo_n_for_l + nr_lo), (0, 0))
                for mO in range(highest_mO):
                    self.add_lo(l_val, (first_helo_n_for_l + nr_lo,) * 2, (mO, mO + 1))

    def to_xml(self) -> ElementTree.Element:
        """Converts the class attributes into an XML structure using ElementTree.

        :return: An ElementTree.Element representing the root of the XML structure.
        """
        spdb = ElementTree.Element("spdb")
        sp = ElementTree.SubElement(spdb, "sp", {k: str(v) for k, v in self.species.items()})
        ElementTree.SubElement(sp, "muffinTin", {k: str(v) for k, v in self.muffin_tin.items()})
        for state in self.atomic_states:
            ElementTree.SubElement(sp, "atomicState", {k: str(v).lower() for k, v in state.items()})

        basis = ElementTree.SubElement(sp, "basis")
        for default in self.basis["default"]:
            ElementTree.SubElement(basis, "default", {k: str(v).lower() for k, v in default.items()})

        for custom in self.basis["custom"]:
            ElementTree.SubElement(basis, "custom", {k: str(v).lower() for k, v in custom.items()})

        for lo in self.basis["lo"]:
            lo_element = ElementTree.SubElement(basis, "lo", l=str(lo["l"]))
            for wf in lo["wf"]:
                ElementTree.SubElement(lo_element, "wf", {k: str(v).lower() for k, v in wf.items()})

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
