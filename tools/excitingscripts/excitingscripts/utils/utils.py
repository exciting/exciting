"""General utils for exciting scripts."""

import re
from os.path import join
from pathlib import Path
from typing import Tuple, List, TypeVar, Union, Dict

import numpy as np
from excitingtools import parse
from scipy.optimize import leastsq, OptimizeResult

# Types for static type checking to support maintaining type consistency. For example, the function
# `sort_lists_by_first_list` should return a tuple of lists with elements belonging to the same type as the elements in
# the lists passed as arguments.
T1 = TypeVar("T1")
T2 = TypeVar("T2")


def sort_lists_by_first_list(
    first_list: List[T1], second_list: List[T2]
) -> Tuple[List[T1], List[T2]]:
    """Sorts two lists, using the first list as reference

    :param first_list: first list to be sorted, used as reference
    :param second_list: second list to be sorted, uses first list as reference
    :return: sorted lists
    """
    first_len = len(first_list)
    second_len = len(second_list)
    assert (
        first_len == second_len
    ), f"Both lists should have the same length, not {first_len} and {second_len}."

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


def get_prettified_scientific_notation(
    number: float, unit: Union[str, None] = None
) -> str:
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
    numbers = re.findall(r"[-+]?\d*\.\d+|\d+", line)
    return [float(num) for num in numbers]


def get_num_atoms(run_dir: str) -> int:
    """Extract the total number of atoms per unit cell from INFO.OUT.

    :param run_dir: directory where exciting runs.
    :return: number of atoms per unit cell.
    """
    # Define the path to the INFO.OUT file
    info_path = join(run_dir, "INFO.OUT")

    # Parsing INFO.OUT using excitingtools
    parsed_info = parse(info_path)

    try:
        return parsed_info["initialization"]["Total number of atoms per unit cell"]
    except KeyError:
        raise ValueError("Number of atoms not found in INFO.OUT")


def get_structure_optimizations_properties(run_dir: str, key: str) -> List[Dict]:
    """Read all lines from the INFO.OUT file, extract property for each optimization step.

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
    """Check the coordinate type is cartesian from input.xml.

    :param run_dir: directory where exciting runs
    :return: coordinate type, either True for "cartesian" or False for "lattice" or other type.
    """
    # Define the path to the INFO.OUT file
    input_path = join(run_dir, "input.xml")

    # Parse the input.xml file
    input_parsed = parse(input_path)

    # Check for cartesian attribute
    return (
        "cartesian" in input_parsed["structure"].keys()
        and input_parsed["structure"]["cartesian"]
    )


def parse_energy_vs_volume(directory: str) -> Tuple[List[float], List[float]]:
    """Read the "energy_vs_volume" file

    :param directory: directory containing the file
    :return: volume and energy values
    """
    volume = []
    energy = []

    data_file = Path(directory) / "energy-vs-volume"

    lines = data_file.read_text().split("\n")

    for line in lines:
        if len(line) != 0:
            vol, ene = line.split()
            volume.append(float(vol))
            energy.append(float(ene))

    return volume, energy


def birch_murnaghan_eos(v: float, p: List[float]) -> float:
    """Calculate energy using the Birch-Murnaghan equation of state.

    :param v: Volume.
    :param p: Parameters for the equation [eq_vol, min_energy, bulk_modulus, bulk_modulus_pressure_deriv].
    :return: Calculated energy.
    """
    eq_vol = p[0]  # Equilibrium volume
    min_energy = p[1]  # Minimum energy
    bulk_modulus = p[2]  # Bulk modulus at equilibrium volume
    bulk_modulus_pressure_deriv = p[3]  # Pressure derivative of bulk modulus

    vol_ratio = (eq_vol / v) ** (2.0 / 3)
    vol_ratio_offset = vol_ratio - 1.0
    energy = min_energy + 9.0 * bulk_modulus * eq_vol / 16.0 * (
        vol_ratio_offset**3 * bulk_modulus_pressure_deriv
        + (6.0 - 4.0 * vol_ratio) * vol_ratio_offset**2
    )

    return energy


def pressure_birch_murnaghan_eos(v: float, p: List[float]) -> float:
    """Calculate pressure using the Birch-Murnaghan equation of state.

    :param v: Volume.
    :param p: Parameters for the equation [eq_vol, min_energy, bulk_modulus, bulk_modulus_pressure_deriv].
    :return: Calculated pressure.
    """
    eq_vol = p[0]  # Equilibrium volume
    bulk_modulus = p[2]  # Bulk modulus at equilibrium volume
    bulk_modulus_pressure_deriv = p[3]  # Pressure derivative of bulk modulus

    vol_ratio = eq_vol / v
    v7 = vol_ratio ** (7.0 / 3)
    v5 = vol_ratio ** (5.0 / 3)
    v2 = vol_ratio ** (2.0 / 3)
    vol_ratio_offset = v2 - 1.0
    pressure_derivative_offset = bulk_modulus_pressure_deriv - 4.0
    pressure = (
        3.0
        * bulk_modulus
        / 2.0
        * (v7 - v5)
        * (1.0 + 3.0 / 4 * pressure_derivative_offset * vol_ratio_offset)
    )
    return pressure


def initial_guess(volumes: List[float], energies: List[float]) -> List[float]:
    """Generate initial guess parameters for the Birch-Murnaghan EOS fit.

    :param volumes: List of volumes.
    :param energies: List of energies.
    :return: Initial guess parameters [eq_vol, min_energy, bulk_modulus, bulk_modulus_pressure_deriv].
    """
    a, b, c = np.polyfit(volumes, energies, 2)  # ax^2 + bx + c
    eq_vol = -b / (2 * a)
    min_energy = a * eq_vol**2 + b * eq_vol + c
    bulk_modulus = 2 * a * eq_vol
    bulk_modulus_pressure_deriv = 2.0
    return [eq_vol, min_energy, bulk_modulus, bulk_modulus_pressure_deriv]


def residuals(p: List[float], e: float, v: float) -> float:
    """Calculate residuals for the least squares fit.

    :param p: Parameters for the Birch-Murnaghan EOS.
    :param e: Energies.
    :param v: Volumes.
    :return: Residuals.
    """
    return e - birch_murnaghan_eos(v, p)


def birch_murnaghan_fit(volumes: List[float], energies: List[float]) -> OptimizeResult:
    """Perform the least squares fit using the Birch-Murnaghan EOS.

    :param volumes: List of volumes.
    :param energies: List of energies.
    :return: Optimized parameters.
    """
    p = initial_guess(volumes, energies)
    return leastsq(residuals, p, args=(energies, volumes))
