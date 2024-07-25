"""Python script for fitting energy-vs-volume curves"""

import os
from argparse import ArgumentParser
from os.path import join
from typing import Tuple, List, Dict, Callable

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ptk
import numpy as np
from excitingscripts.utils.utils import (
    sort_lists_by_first_list,
    parse_energy_vs_volume,
    birch_murnaghan_eos,
    pressure_birch_murnaghan_eos,
    residuals,
    birch_murnaghan_fit,
)
from scipy.constants import physical_constants

if matplotlib.__version__.split(".")[0] == "2":
    matplotlib.style.use("classic")

symmetry_factors = {0: 1, 1: 1, 2: 2, 3: 4}
symmetry_labels = {0: "     ", 1: "(sc) ", 2: "(bcc)", 3: "(fcc)"}

# Unit conversion factor setup
bohr_radius = physical_constants["Bohr radius"][0]  # in meters
joule2hartree = physical_constants["hartree-joule relationship"][0]  # in joules

unitconv = joule2hartree / (bohr_radius**3) * 10**-9

fmt = "%10.4f"
amt = "%10.4f"
bmt = "%8.3f"
lmt = "%10.2f"


def pressure_finite_difference(
    volumes: List[float], energies: List[float]
) -> Tuple[List[float], List[float]]:
    """Calculate pressure using the Finite Difference method

    :param volumes: List of volumes.
    :param energies: List of energies.
    :return: Calculated pressure.
    """
    interpolated_vol = []
    pressure_differences = []
    for i in range(len(volumes) - 1):
        pressure_difference = (energies[i + 1] - energies[i]) / (
            volumes[i] - volumes[i + 1]
        )
        pressure_differences.append(pressure_difference * unitconv * 10)
        interpolated_vol.append(volumes[i] + (volumes[i + 1] - volumes[i]) / 2)
    return interpolated_vol, pressure_differences


def fit_pressure_vs_volume(
    volumes: List[float], energies: List[float]
) -> Tuple[Callable, List[float], List[float], List[float], Dict, float]:
    """
    Fit the energy-vs-volume data using the Birch-Murnaghan Equation of State.

    :param volumes: List of volumes.
    :param energies: List of energies.

    :return:Birch fit, bulk modulus, bulk_modulus_pressure_deriv, minima, lattice constant, chi value
    """
    p = birch_murnaghan_fit(volumes, energies)[0]
    eq_vol, min_energy, bulk_modulus, bulk_modulus_pressure_deriv = p
    bulk_modulus = bulk_modulus * unitconv

    lattice_const = {i: [] for i in [1, 2, 3]}

    for isym in [1, 2, 3]:
        lattice_const[isym].append(np.cbrt(symmetry_factors[isym] * eq_vol))

    chi = 0
    energy_bm_eos = []
    fitted_chi = 0
    for i in range(len(volumes)):
        chi = chi + residuals(p, energies[i], volumes[i]) ** 2
        fitted_chi = (
            fitted_chi
            + (
                energies[i]
                - birch_murnaghan_eos(
                    volumes[i],
                    [
                        eq_vol,
                        min_energy,
                        bulk_modulus / unitconv,
                        bulk_modulus_pressure_deriv,
                    ],
                )
            )
            ** 2
        )
        energy_bm_eos.append(
            birch_murnaghan_eos(
                volumes[i],
                [
                    eq_vol,
                    min_energy,
                    bulk_modulus / unitconv,
                    bulk_modulus_pressure_deriv,
                ],
            )
        )

    chi = np.sqrt(chi) / len(volumes)

    def curv(xvol):
        fitted_pressures = []
        for x in xvol:
            fitted_pressures.append(
                pressure_birch_murnaghan_eos(
                    x,
                    [
                        eq_vol,
                        min_energy,
                        bulk_modulus * 10,
                        bulk_modulus_pressure_deriv,
                    ],
                )
            )
        return fitted_pressures

    return (
        curv,
        [bulk_modulus],
        [bulk_modulus_pressure_deriv],
        [eq_vol],
        lattice_const,
        chi,
    )


def print_info_to_stdout(
    dmin: List[float],
    lattice_const: Dict,
    bulk_modulus: List[float],
    bulk_modulus_pressure_deriv: List[float],
    chi: float,
) -> None:
    """Print information onto the screen

    :param dmin: minima of the fit.
    :param lattice_const: lattice constant for the given lattice symmetry code
    :param bulk_modulus: bulk modulus
    :param bulk_modulus_pressure_deriv: derivative of bulk modulus with respect to pressure
    :param chi: chi-squared value indicating the goodness of fit
    """
    if len(dmin) > 1:
        print("##############################################\n")
        print("WARNING: Multiple minima are found!\n")
        print("##############################################\n")

    if len(dmin) == 0:
        print("##############################################\n")
        print("WARNING: No minimum in the given xrange!\n")
        print("##############################################\n")

    print(
        "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    )
    print(
        "     V0        B0         Bp        a-sc       a-bcc      a-fcc     log(chi)"
    )
    for i, eq_vol in enumerate(dmin):
        print(
            fmt % eq_vol,
            bmt % (bulk_modulus[i]),
            bmt % (bulk_modulus_pressure_deriv[i]),
            amt % (lattice_const[1][i]),
            amt % (lattice_const[2][i]),
            amt % (lattice_const[3][i]),
            lmt % (np.log10(chi)),
        )
    print(
        "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    )


def plot_pressure_vs_volume(
    volumes: List[float],
    volumes_fd,
    pressure_fd,
    curv: np.poly1d,
    dmin: np.ndarray,
    output_dir: str,
):
    """
    Plot the pressure-vs-volume data and save the plot.

    :param volumes : List of volumes.
    :param volumes_fd : List of volumes at which pressure is calculated using finite difference.
    :param pressure_fd: List of pressures calculated using finite difference.
    :param curv : Fitted curve.
    :param dmin : Minima of the fit.
    :param output_dir : Directory to save the plot.
    """
    plt.rcParams.update(
        {
            "figure.figsize": (10, 7.5),
            "ytick.minor.size": 6,
            "xtick.major.pad": 8,
            "ytick.major.pad": 4,
            "patch.linewidth": 2.0,
            "axes.linewidth": 2.0,
            "lines.linewidth": 3.5,
            "lines.markersize": 12.0,
            "axes.formatter.limits": (-4, 6),
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "axes.labelsize": 20,
            "axes.labelcolor": "black",
            "axes.axisbelow": "True",
            "legend.fontsize": 22,
            "axes.titlesize": 24,
            "axes.titlepad": 20,
        }
    )
    # Create figure and axis
    fig, ax = plt.subplots()

    yfmt = ptk.ScalarFormatter(useOffset=True, useMathText=True)

    # Set labels
    ax.set_xlabel("Volume [Bohr\u00b3]")
    ax.set_ylabel(r"Pressure [kbar]")

    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(2)

    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.grid(True)

    # Plotting the pressure vs volume
    xvol = np.linspace(min(volumes), max(volumes), 100)
    plt.plot(xvol, curv(xvol), "b-", label="birch-murnaghan fit")
    plt.plot(volumes_fd, pressure_fd, "go", label="finite differences")

    plt.legend(loc=9, borderaxespad=0.8, numpoints=1)

    ymax = max(max(curv(xvol)), max(pressure_fd))
    ymin = min(min(curv(xvol)), min(pressure_fd))
    dxx = abs(max(xvol) - min(xvol)) / 18
    dyy = abs(ymax - ymin) / 18
    ax.yaxis.set_major_formatter(yfmt)
    ax.set_xlim(min(xvol) - dxx, max(xvol) + dxx)
    ax.set_ylim(ymin - dyy, ymax + dyy)

    ax.xaxis.set_major_locator(ptk.MaxNLocator(7))

    plt.tight_layout()
    plt.savefig(
        join(output_dir, "PLOT.png"), orientation="portrait", format="png", dpi=300
    )

    plt.close()


def main() -> None:
    parser = ArgumentParser(
        description="Python script for fitting energy-vs-volume curves."
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

    args = parser.parse_args()
    directory = args.directory[0]

    # Read data from energy-vs-volume file
    volumes, energies = parse_energy_vs_volume(directory)
    volumes, energies = sort_lists_by_first_list(volumes, energies)

    # Fit the data to the given order of polynomial
    curv, bulk_modulus, bulk_modulus_pressure_deriv, dmin, lattice_const, chi = (
        fit_pressure_vs_volume(volumes, energies)
    )

    # Print results
    print_info_to_stdout(
        dmin, lattice_const, bulk_modulus, bulk_modulus_pressure_deriv, chi
    )

    # Get pressure calculated from finite differences
    volumes_fd, pressure_fd = pressure_finite_difference(volumes, energies)

    # Plot results
    plot_pressure_vs_volume(volumes, volumes_fd, pressure_fd, curv, dmin, directory)


if __name__ == "__main__":
    main()
