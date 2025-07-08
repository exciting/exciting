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


def vinet_eos(
    v: float,
    eq_vol: float,
    min_energy: float,
    bulk_modulus: float,
    bulk_modulus_pressure_deriv: float,
) -> float:
    """
    Vinet equation of state.

    :param v: Volume
    :param eq_vol: Equilibrium volume
    :param min_energy: Minimum energy
    :param bulk_modulus: Bulk modulus at equilibrium volume
    :param bulk_modulus_pressure_deriv: Pressure derivative of bulk modulus
    :return: Energy at volume v
    """
    vol_ratio = (v / eq_vol) ** (1.0 / 3)
    gruneisen_param = 3.0 / 2.0 * (bulk_modulus_pressure_deriv - 1.0)
    energy = min_energy - 9 * eq_vol * bulk_modulus * np.exp(
        gruneisen_param * (1 - vol_ratio)
    ) * (1.0 / gruneisen_param**2 - 1.0 / gruneisen_param + vol_ratio / gruneisen_param)
    return energy


def fit_energy_vs_volume(
    volumes: List[float], energies: List[float]
) -> Tuple[Callable, List[float], List[float], List[float], Dict, float]:
    """
    Fit the energy-vs-volume data using the Vinet Equation of State.

    :param volumes: List of volumes.
    :param energies: List of energies.

    :return:Vinet fit, bulk modulus, bulk_modulus_pressure_deriv, minima, lattice constant, chi value
    """
    p = birch_murnaghan_fit(volumes, energies)[0]
    eq_vol, min_energy, bulk_modulus, bulk_modulus_pressure_deriv = p
    bulk_modulus = bulk_modulus * unitconv

    lattice_const = {i: [] for i in [1, 2, 3]}

    for isym in [1, 2, 3]:
        lattice_const[isym].append(np.cbrt(symmetry_factors[isym] * eq_vol))

    chi = 0
    energy_vinet_eos = []
    fitted_chi = 0
    for i in range(len(volumes)):
        chi = chi + residuals(p, energies[i], volumes[i]) ** 2
        fitted_chi = (
            fitted_chi
            + (
                energies[i]
                - vinet_eos(
                    volumes[i],
                    eq_vol,
                    min_energy,
                    bulk_modulus / unitconv,
                    bulk_modulus_pressure_deriv,
                )
            )
            ** 2
        )
        energy_vinet_eos.append(
            vinet_eos(
                volumes[i],
                eq_vol,
                min_energy,
                bulk_modulus / unitconv,
                bulk_modulus_pressure_deriv,
            )
        )

    chi = np.sqrt(chi) / len(volumes)

    def curv(xvol):
        fitted_energies = []
        for x in xvol:
            fitted_energies.append(
                vinet_eos(
                    x,
                    eq_vol,
                    min_energy,
                    bulk_modulus / unitconv,
                    bulk_modulus_pressure_deriv,
                )
            )
        return fitted_energies

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


def plot_energy_vs_volume(
    volumes: List[float],
    energies: List[float],
    curv: np.poly1d,
    dmin: np.ndarray,
    output_dir: str,
):
    """
    Plot the energy-vs-volume data and save the plot.

    :param volumes : List of volumes.
    :param energies : List of energies.
    :param curv : Polynomial fit curve.
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
    ax.set_xlabel("Volume [Bohr³]")
    ax.set_ylabel("Energy [Ha]")

    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(2)

    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.grid(True)

    xvol = np.linspace(min(volumes), max(volumes), 100)
    plt.plot(xvol, curv(xvol), "b-", label=f"vinet fit")
    plt.plot(volumes, energies, "go", label="calculated")
    plt.plot(dmin, curv(dmin), "ro")
    plt.legend(loc=9, borderaxespad=0.8, numpoints=1)

    ymax = max(max(curv(xvol)), max(energies))
    ymin = min(min(curv(xvol)), min(energies))
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
        fit_energy_vs_volume(volumes, energies)
    )

    # Print results
    print_info_to_stdout(
        dmin, lattice_const, bulk_modulus, bulk_modulus_pressure_deriv, chi
    )

    # Plot results
    plot_energy_vs_volume(volumes, energies, curv, dmin, directory)


if __name__ == "__main__":
    main()
