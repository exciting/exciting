"""Python script for fitting energy-vs-volume curves"""

import os
from argparse import ArgumentParser
from os.path import join
from typing import Tuple, List, Dict, Callable

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ptk
import numpy as np
from excitingscripts.plot.pbirch import print_info_to_stdout
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


def bulk_modulus_finite_difference(
    volumes: List[float], energies: List[float]
) -> Tuple[List[float], List[float]]:
    """Calculate bulk modulus using the Finite Difference method

    :param volumes: List of volumes.
    :param energies: List of energies.
    :return: Calculated bulk modulus.
    """
    interpolated_vol_pressure = []
    pressure_differences = []
    for i in range(len(volumes) - 1):
        pressure_difference = (energies[i + 1] - energies[i]) / (
            volumes[i] - volumes[i + 1]
        )
        pressure_differences.append(pressure_difference * unitconv * 10)
        interpolated_vol_pressure.append(volumes[i] + (volumes[i + 1] - volumes[i]) / 2)

    interpolated_vol_bulk_modulus = []
    bulk_modulus_differences = []
    for i in range(len(interpolated_vol_pressure) - 1):
        bulk_modulus_gradient = (
            pressure_differences[i + 1] - pressure_differences[i]
        ) / (interpolated_vol_pressure[i] - interpolated_vol_pressure[i + 1])
        midpoint_volume = (
            interpolated_vol_pressure[i]
            + (interpolated_vol_pressure[i + 1] - interpolated_vol_pressure[i]) / 2
        )
        bulk_modulus_differences.append(bulk_modulus_gradient * midpoint_volume)
        interpolated_vol_bulk_modulus.append(midpoint_volume)

    return interpolated_vol_bulk_modulus, bulk_modulus_differences


def fit_pressure_vs_volume(
    volumes: List[float], energies: List[float]
) -> Tuple[Callable, List[float], List[float], List[float], Dict, float]:
    """Fit the energy-vs-volume data using the Birch-Murnaghan Equation of State.

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
        chi += residuals(p, energies[i], volumes[i]) ** 2
        fitted_chi += (
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
        ) ** 2
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
        volume_midpoints = []
        calculated_bulk_modulus = []
        for i in range(len(fitted_pressures) - 1):
            pressure_gradient = (fitted_pressures[i + 1] - fitted_pressures[i]) / (
                xvol[i] - xvol[i + 1]
            )
            midpoint_volume = xvol[i] + (xvol[i + 1] - xvol[i]) / 2
            calculated_bulk_modulus.append(pressure_gradient * midpoint_volume)
            volume_midpoints.append(midpoint_volume)
        return volume_midpoints, calculated_bulk_modulus

    return (
        curv,
        [bulk_modulus],
        [bulk_modulus_pressure_deriv],
        [eq_vol],
        lattice_const,
        chi,
    )

def findindex(x: float, y: List[float], dymax=1e10) -> int:
    """Finds the index of given value in the list upto specified tolerance

    :param x: value for which index is found
    :param y: list in which the value is searched for
    :param dymax: tolerance for which value can be found
    :return : index of the value in the list
    """
    k = 0
    for j in range(len(y)):
        if abs(x - y[j]) < dymax:
            k = j
            dymax = abs(x - y[j])
    return k


def plot_bulk_modulus_vs_volume(
    volumes: List[float],
    volumes_fd,
    bulk_modulus_fd,
    curv: np.poly1d,
    dmin: np.ndarray,
    output_dir: str,
):
    """
    Plot the bulk modulus vs volume data and save the plot.

    :param volumes : List of volumes.
    :param volumes_fd : List of volumes at which bulk modulus is calculated using finite difference.
    :param bulk_modulus_fd: List of bulk modulus calculated using finite difference.
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
            "lines.linewidth": 2.5,
            "lines.markersize": 10.0,
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
    ax.set_ylabel(r"Bulk modulus [kbar]")

    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(2)

    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.grid(True)

    # Plotting the pressure vs volume
    xvol = np.linspace(min(volumes), max(volumes), 100)
    bvol, bulk_modulus = curv(xvol)

    # Calculating and printing chi-value for Bulk Modulus
    bchi = 0
    for i in range(len(volumes_fd)):
        k = findindex(volumes_fd[i], bvol)
        bchi = bchi + (bulk_modulus[k] - bulk_modulus_fd[i]) ** 2

    bchi = np.sqrt(bchi / len(volumes_fd)) / 10.0

    print("\n     DB0 =", bmt % bchi, "\n")

    plt.plot(bvol, bulk_modulus, "b-", label="birch-murnaghan fit")
    plt.plot(volumes_fd, bulk_modulus_fd, "go", label="finite differences")

    plt.legend(loc=9, borderaxespad=0.8, numpoints=1)

    ymax = max(bulk_modulus)
    ymin = min(bulk_modulus)
    dxx = abs(max(bvol) - min(bvol)) / 18
    dyy = abs(ymax - ymin) / 18
    ax.yaxis.set_major_formatter(yfmt)
    ax.set_xlim(min(bvol) - dxx, max(bvol) + dxx)
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
    volumes_fd, pressure_fd = bulk_modulus_finite_difference(volumes, energies)

    # Plot results
    plot_bulk_modulus_vs_volume(volumes, volumes_fd, pressure_fd, curv, dmin, directory)


if __name__ == "__main__":
    main()
