"""Check-fit implementation."""

import json
import math
import warnings
from argparse import ArgumentParser
from pathlib import Path
from typing import Optional, List, Callable

import numpy as np
from scipy.constants import physical_constants

from excitingscripts.utils.utils import sort_lists_by_first_list

warnings.filterwarnings("error")
invcm2hz = physical_constants["hertz-inverse meter relationship"][0] / 100

try:
    # Check if np.RankWarning exists directly
    numpy_RankWarning = np.RankWarning
except AttributeError:
    # If not, fall back to np.exceptions.RankWarning
    numpy_RankWarning = np.exceptions.RankWarning

def arg_parser(quantity: str) -> ArgumentParser:
    """Get the arg parser for checkfit scripts.

    :param quantity: what derivatives to check fit
    :return: the argparser
    """
    parser = ArgumentParser(
        description=f"Python script for extracting the derivatives at zero displacement of the "
        f"quantity-vs-displacement curves for quantity: {quantity}."
    )

    parser.add_argument(
        "maximum_displacement_fit", type=float, help="Maximum displacement in the unit cell used for the fit."
    )

    parser.add_argument("order_of_derivative", type=int, help="Order of derivative for the fit.")

    parser.add_argument("atomic_mass", type=float, help="Atomic mass in atomic-mass-units.")

    return parser


def fit(order: int, x: List[float], y: List[float], order_of_derivative: int) -> Optional[float]:
    """Fit data to a polynomial.

    :param order: order of the fit
    :param x: x data
    :param y: y data
    :param order_of_derivative: determines which coefficient is taken from the fit
    :return: the chosen coefficient, if the fit is poorly conditioned "None"
    """
    try:
        fit_coefficient = np.polyfit(x, y, order)
    except numpy_RankWarning:
        return None
    return math.factorial(order_of_derivative) * fit_coefficient[order - order_of_derivative]


def print_info_to_stdout(
    max_displacement: float, frequencies: List[float], order_of_derivative: int, n_max: int
) -> None:
    """Print some information to the terminal.

    :param max_displacement: maximum chosen displacement for the fit
    :param frequencies: the computed frequencies
    :param order_of_derivative: fit order of interest
    :param n_max: number of chosen displacement values
    """
    text = (
        f"\n Fit data\n\n"
        f" Maximum value of the displacement: {max_displacement}\n"
        f" Number of displacement values used: {n_max}\n"
        f" Fit results for the derivative of order: {order_of_derivative}\n\n"
    )

    not_none_frequencies = filter(lambda x: x is not None, frequencies)
    hertz_text = ""
    for j, frequency in enumerate(not_none_frequencies):
        text += f" Polynomial of order  {order_of_derivative + j}  ==>   {frequency:.2f} [cm-1]\n"
        freq_in_thz = frequency / invcm2hz / 10**12
        hertz_text += f" Polynomial of order  {order_of_derivative + j}  ==>   {freq_in_thz:.4f} [THz]\n"
    print(text + "\n" + hertz_text + "\n")


def get_unit_conversion_factor(lattice_param_exp: int) -> float:
    """Get the unit conversion factor.

    :param lattice_param_exp: exponent for the lattice parameter in the unit conversion.
    :return: unit conversion factor
    """
    info_diamond_phonon = Path("INFO-diamond-phonon.json")
    if not info_diamond_phonon.exists():
        raise FileNotFoundError(f"file {info_diamond_phonon.name} not found")

    with open(info_diamond_phonon) as fid:
        info = json.load(fid)
    lattice_parameter = info["Equilibrium lattice parameter (alat) in a.u."]

    amu = physical_constants["atomic mass constant"][0]
    bohr_radius = physical_constants["Bohr radius"][0]
    joule2hartree = physical_constants["hartree-joule relationship"][0]

    return (invcm2hz / bohr_radius) ** 2 * joule2hartree / (amu * lattice_parameter**lattice_param_exp)


def quantity_specific_checkfit(
    quantity: str, factor: float, lattice_param_exp: int
) -> Callable[[float, int, float], None]:
    """Specific checkfit implementation for either energy or force.

    :param quantity: name of the quantity key in the result file
    :param factor: factor for the unit conversion
    :param lattice_param_exp: exponential for the lattice parameter in the unit conversion
    return: checkfit function
    """

    def checkfit(max_displacement_fit: float, order_of_derivative: int, atomic_mass: float) -> None:
        """Fit the computed values to some polynomial function.

        Note: Assumes that we have a symmetric displacement list, each value and its negative value should be in it

        :param max_displacement_fit: the maximum displacement taken for the fit
        :param order_of_derivative: determines which coefficient is taken from the fit
        :param atomic_mass: the mass in atomic units from the atom-species in the unit cell
        """
        result_file = Path("phonon_results.json")
        if not result_file.exists():
            raise FileNotFoundError(f"file {result_file.name} not found")

        assert order_of_derivative >= 0, "Order of derivative must be positive"

        full_conversion_factor = get_unit_conversion_factor(lattice_param_exp) * factor / atomic_mass

        with open(result_file) as fid:
            phonon_data: dict = json.load(fid)["results"]

        quantity_values = []
        displacements = []
        for displacement_string, result in phonon_data.items():
            displacement = float(displacement_string)
            if abs(displacement) <= max_displacement_fit:
                quantity_values.append(result[quantity])
                displacements.append(displacement)

        displacements, quantity_values = sort_lists_by_first_list(displacements, quantity_values)

        n_max = len(displacements)
        border = max(order_of_derivative, 1)
        number_of_fits = 6
        results = []

        while len(displacements) > border:
            frequencies = []
            max_displacement = max(displacements)
            for order in range(order_of_derivative, order_of_derivative + number_of_fits):
                fit_coefficient = fit(order, displacements, quantity_values, order_of_derivative)
                frequency = None
                if fit_coefficient is not None:
                    frequency = np.sqrt(abs(fit_coefficient) * full_conversion_factor) / (2 * np.pi)
                frequencies.append(frequency)

            if len(displacements) == n_max:
                print_info_to_stdout(max_displacement, frequencies, order_of_derivative, n_max)
            results.append({"max_displacement": max_displacement, "frequencies": frequencies})

            displacements = displacements[1:-1]
            quantity_values = quantity_values[1:-1]

        full_results = {"order_of_derivative": order_of_derivative, "fits": results}
        with open(f"checkfit_{quantity}_results.json", "w") as fid:
            json.dump(full_results, fid, indent=4)

    return checkfit


if __name__ == "__main__":
    print("This is no script to call directly. Please call a quantity-specific version!")
