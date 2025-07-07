"""Python script for extracting derivatives at zero strain of energy-vs-strain curves."""

import json
import os
from argparse import ArgumentParser
from os.path import join
from pathlib import Path
from typing import Tuple, List, Dict

from excitingscripts.checkfit.checkfit import fit
from excitingscripts.utils.utils import sort_lists_by_first_list
from scipy.constants import physical_constants


def parse_info_elastic_constants(directory=os.getcwd()) -> Dict:
    """Gives the necessary info to be printed

    :param directory: path to the directory containing the file "INFO-elastic-constants"
    :return : dictionary containing the information
    """
    info_file = Path(directory) / "INFO-elastic-constants"

    info = {
        "Maximum Lagrangian strain": None,
        "Number of strain values": None,
        "Volume of equilibrium unit cell": None,
        "Deformation code": None,
        "Deformation label": None,
    }

    lines = info_file.read_text().split("\n")

    for line in lines:
        if "Maximum Lagrangian strain" in line:
            info["Maximum Lagrangian strain"] = float(line.split("=")[-1])
        elif "Number of strain values" in line:
            info["Number of strain values"] = int(line.split("=")[-1])
        elif "Volume of equilibrium unit cell" in line:
            info["Volume of equilibrium unit cell"] = float(
                line.split("=")[-1].split()[0]
            )
        elif "Deformation code" in line:
            info["Deformation code"] = line.split("=")[-1].strip()
        elif "Deformation label" in line:
            info["Deformation label"] = line.split("=")[-1].strip()

    return info


def parse_energy_vs_strain(
    directory: str, max_strain: float
) -> Tuple[List[float], List[float]]:
    """Read the "energy_vs_strain" file

    :param directory: directory containing the file
    :param max_strain: value of maximum strain
    :return: strain and energy values
    """
    strain = []
    energy = []

    data_file = join(directory, "energy-vs-strain")

    if not os.path.exists(data_file):
        raise FileNotFoundError("energy-vs-strain not found")

    with open(data_file, "r") as input_energy:
        for line in input_energy:
            eta, ene = map(float, line.split())
            if abs(eta) <= max_strain:
                strain.append(eta)
                energy.append(ene)
    return strain, energy


def print_info_to_stdout(
    info: dict,
    max_strain: float,
    derivatives: List[Dict[str, float]],
    order_of_derivative: int,
    n_max: int,
) -> None:
    """Print some information to the terminal.

    :param info:
    :param max_strain: maximum chosen strain for the fit
    :param derivatives: the computed derivatives
    :param order_of_derivative: fit order of interest
    :param n_max: number of chosen strain values
    """
    print("\n###########################################\n")
    print("Fit data-----------------------------------\n")
    print(f"Deformation code             ==> {info['Deformation code']}")
    print(f"Deformation label            ==> {info['Deformation label']}")
    print(f"Maximum value of the strain  ==> {max_strain:.8f}")
    print(f"Number of strain values used ==> {n_max}\n")
    print(f"Fit results for the derivative of order   {order_of_derivative}\n")
    for derivative in derivatives:
        print(
            f"Polynomial of order  {derivative['order']} ==>   {derivative['value']:.2f} [GPa]"
        )
    print("\n###########################################\n")


def save_derivatives_to_json(
    directory: str, order_of_derivative: int, results: List[Dict[str, float]]
) -> None:
    """Saving the data in fitted derivative data in a JSON format

    :param directory: path to save the file
    :param order_of_derivative: order of derivative for the fit
    :param results: results of the fits
    """
    full_results = {"order_of_derivative": order_of_derivative, "fits": results}

    with open(join(directory, "checkfit_energy_results.json"), "w") as f:
        json.dump(full_results, f, indent=4)


def main() -> None:
    parser = ArgumentParser(
        description="Python script for extracting derivatives at zero strain of energy-vs-strain curves."
    )

    parser.add_argument(
        "--root-directory",
        "-r",
        type=str,
        default=[os.getcwd()],
        nargs=1,
        dest="directory",
        help="Directory containing the necessary files",
    )

    parser.add_argument(
        "--maximum-strain",
        "--s-max",
        type=float,
        required=True,
        nargs=1,
        dest="max_strain",
        help="Maximum strain for the fit",
    )

    parser.add_argument(
        "--order-of-derivative",
        "-o",
        type=int,
        required=True,
        nargs=1,
        dest="order_of_derivative",
        help="The order of derivative",
    )

    args = parser.parse_args()
    directory = args.directory[0]
    max_strain = args.max_strain[0]
    order_of_derivative = args.order_of_derivative[0]

    # Parse INFO-elastic-constants file
    info = parse_info_elastic_constants(directory)

    # Read and filter energy-vs-strain data
    strains, energy = parse_energy_vs_strain(directory, max_strain)

    # Sort data
    strains, energy = sort_lists_by_first_list(strains, energy)

    # Prepare polynomial orders
    orderlist = [order_of_derivative + i for i in range(6)]

    # Unit conversion factor setup
    bohr_radius = physical_constants["Bohr radius"][0]  # in meters
    joule2hartree = physical_constants["hartree-joule relationship"][0]  # in joules
    volume = info["Volume of equilibrium unit cell"]

    unitconv = joule2hartree / (bohr_radius**3 * volume) * 10**-9

    # Fit data and calculate derivatives
    n_max = len(strains)
    results = []
    while len(strains) > order_of_derivative and len(strains) > 1:
        derivatives = [
            fit(order, strains, energy, order_of_derivative) for order in orderlist
        ]
        derivatives = [d * unitconv if d is not None else None for d in derivatives]

        strain = max(strains)
        results.append(
            {
                "max_strain": strain,
                "derivatives": [
                    {"order": order_of_derivative + i, "value": deriv}
                    for i, deriv in enumerate(derivatives)
                    if deriv is not None
                ],
            }
        )

        strains = strains[1:-1]
        energy = energy[1:-1]

    # Print results
    print_info_to_stdout(
        info, max_strain, results[0]["derivatives"], order_of_derivative, n_max
    )

    # Save derivatives to JSON file
    save_derivatives_to_json(directory, order_of_derivative, results)


if __name__ == "__main__":
    main()
