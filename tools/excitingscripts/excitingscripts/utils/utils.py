"""General utils for exciting scripts."""
import numpy as np
import re
import os
from os.path import join, exists
from pathlib import Path
from typing import Tuple, List, TypeVar, Union, Dict
from excitingtools import parse

# Types for static type checking to support maintaining type consistency. For example, the function
# `sort_lists_by_first_list` should return a tuple of lists with elements belonging to the same type as the elements in
# the lists passed as arguments.
T1 = TypeVar("T1")
T2 = TypeVar("T2")


def sort_lists_by_first_list(first_list: List[T1], second_list: List[T2]) -> Tuple[List[T1], List[T2]]:
    """ Sorts two lists, using the first list as reference

    :param first_list: first list to be sorted, used as reference
    :param second_list: second list to be sorted, uses first list as reference
    :return: sorted lists
    """
    first_len = len(first_list)
    second_len = len(second_list)
    assert first_len == second_len, f"Both lists should have the same length, not {first_len} and {second_len}."

    sorted_indices = np.argsort(first_list)
    sorted_first_list = [first_list[x] for x in sorted_indices]
    sorted_second_list = [second_list[x] for x in sorted_indices]

    return sorted_first_list, sorted_second_list


def get_decimal_decomposition(number: float) -> Tuple[float, int]:
    """Decompose the number into mantissa and exponent.

    :param number: input number
    :return: tuple with shifted number (only one leading digit before the decimal point) and exponent
    """
    log_num = np.log10(abs(number))
    exponent = int(log_num)
    shifted_number = 10 ** (log_num - exponent)
    while shifted_number < 1:
        shifted_number *= 10
        exponent -= 1
    return shifted_number, exponent


def get_prettified_scientific_notation(number: float, unit: Union[str, None] = None) -> str:
    """Decompose the number into mantissa and exponent and produce formatted string.

    :param number: input number
    :param unit: unit of the number
    :return: prettified string representation
    """
    shifted_number, exponent = get_decimal_decomposition(number)

    exponent_string = f"$10^{{{exponent}}}$"
    if exponent == 0:
        exponent_string = ""
    elif exponent == 1:
        exponent_string = r"$10\,$"

    sign = "+" if number >= 0 else "\u2013"
    representation = rf"${sign}{shifted_number:6.4f}\cdot${exponent_string}"
    if unit is None:
        return representation
    return representation + f"[{unit}]"


def extract_values_from_line(line: str) -> List[float]:
    """Extract all numbers from a given line using regular expressions.

    :param line: input string from which to extract numbers.
    :return: list of values found in the input string.
    """
    numbers = re.findall(r'[-+]?\d*\.\d+|\d+', line)
    return [float(num) for num in numbers]


def get_num_atoms(run_dir: str) -> int:
    """ Extract the total number of atoms per unit cell from INFO.OUT.

    :param run_dir: directory where exciting runs.
    :return: number of atoms per unit cell.
    """
    # Define the path to the INFO.OUT file
    info_path = join(run_dir, "INFO.OUT")

    # Parsing INFO.OUT using excitingtools
    parsed_info = parse(info_path)

    try:
        return parsed_info['initialization']['Total number of atoms per unit cell']
    except KeyError:
        raise ValueError("Number of atoms not found in INFO.OUT")


def get_structure_optimizations_properties(run_dir: str, key: str) -> List[Dict]:
    """ Read all lines from the INFO.OUT file, extract property for each optimization step.

    :param run_dir: directory where exciting runs.
    :param key: property name which is parsed for each optimization step. Available ones are:
                "Maximum force",
                "Center of mass",
                "Total torque",
                "Number of total scf iterations",
                "Total atomic forces",
                "Total energy",
                "Atomic positions"

    :return: list of dictionaries containing properties.
    """
    # Define the path to the INFO.OUT file
    info_path = join(run_dir, "INFO.OUT")

    data = []

    # Parsing using excitingtools
    parsed_info = parse(info_path)

    for i in parsed_info["str_opt"].keys():
        if key in parsed_info["str_opt"][i].keys():
            data.append(parsed_info["str_opt"][i][key])
        else:
            raise ValueError(f"{key} doesn't exist in INFO.OUT")

    return data


def is_coordinate_cartesian(run_dir: str) -> str:
    """ Check the coordinate type is cartesian from input.xml.

    :param run_dir: directory where exciting runs
    :return: coordinate type, either True for "cartesian" or False for "lattice" or other type.
    """
    # Define the path to the INFO.OUT file
    input_path = join(run_dir, "input.xml")

    # Parse the input.xml file
    input_parsed = parse(input_path)

    # Check for cartesian attribute
    return 'cartesian' in input_parsed['structure'].keys() and input_parsed['structure']['cartesian']
