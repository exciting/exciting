"""Extract the trajectory of the optimization."""

import copy
import re
from pathlib import Path
from typing import Dict, List, Union

import numpy as np

from excitingtools import ExcitingStructure
from excitingtools.exciting_dict_parsers.groundstate_parser import parse_info_out


def parse_initialization_structure(initialization_d: Dict) -> ExcitingStructure:
    """Parse the initialization structure.

    :param initialization_d: dictionary containing the initialization structure, parsed from the INFO.OUT file
    :return: the initialization structure
    """
    lattice = np.reshape(initialization_d["Lattice vectors (cartesian)"], (3, 3), "F")

    pattern = re.compile(r"^Species (\d+)$")
    matches = [k for k in initialization_d if pattern.match(k)]

    my_atoms = []
    for spec_key in matches:
        spec_info = initialization_d[spec_key]
        positions = spec_info["Atomic positions"]
        for str_posi in positions.values():
            my_posi = [float(x) for x in str_posi.split()]
            my_atoms.append({"species": spec_info["Species symbol"], "position": my_posi})

    return ExcitingStructure(my_atoms, lattice)


def parse_optimization_trajectory(info_out_file: Union[str, Path]) -> List[ExcitingStructure]:
    """Parse the optimization trajectory.

    Read the INFO.OUT to find the optimization trajectory.

    :param info_out_file: path to the INFO.OUT file
    :return: list of structures
    """
    info_dict = parse_info_out(str(info_out_file))

    reference_structure = parse_initialization_structure(info_dict["initialization"])
    trajectory = [reference_structure]

    str_opt_info = info_dict.get("str_opt")
    if str_opt_info is None or len(str_opt_info) == 0:
        return trajectory

    # skip the first entry as it is the same as the initial structure
    for i in range(1, len(str_opt_info)):
        step_struct = copy.deepcopy(reference_structure)
        step_struct.positions = list(str_opt_info[i]["Atomic positions"].values())
        trajectory.append(step_struct)

    return trajectory
