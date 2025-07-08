"""Python visualization tool for relaxed coordinates of atoms in the unit cell."""

import os
from argparse import ArgumentParser
from os.path import join
from typing import List, Dict

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ptk
import numpy as np
from excitingscripts.utils.utils import (
    get_structure_optimizations_properties,
    get_num_atoms,
    is_coordinate_cartesian,
)

# Use classic style if Matplotlib version is 2
if matplotlib.__version__.split(".")[0] == "2":
    matplotlib.style.use("classic")


def get_optimized_relative_coordinates(
    atom1: int, atom2: int, root_directory=os.getcwd()
) -> List[float]:
    """Get the optimized relative coordinates between two atoms.

    :param atom1: first atom number.
    :param atom2: second atom number.
    :param root_directory: root directory.
    :return: List of relative coordinates.
    """
    # Counting the number of rundir-xx
    listdir = os.listdir(root_directory)
    displ_points = len(listdir)

    for directory in listdir:
        if not directory.startswith("rundir"):
            displ_points -= 1

    rundir_infty = f"{root_directory}/rundir-oo"
    if os.path.exists(rundir_infty):
        displ_points -= 1

    if displ_points == 0:
        raise ValueError("Directory doesn't have any rundir")

    optimized_geometries = []
    max_strains = []

    for i in range(displ_points + 1):
        if not os.path.exists(rundir_infty) and i == 0:
            continue

        # Specify the directory
        rundir_i = join(root_directory, f"rundir-{i}")

        # Get the optimized coordinates
        # Extract data from the file
        optimized_position_data: List[Dict] = get_structure_optimizations_properties(
            rundir_i, "Atomic positions"
        )[-1]

        # List to store the relative coordinates

        pos1 = optimized_position_data[atom1]
        pos2 = optimized_position_data[atom2]

        # Find the relative position between them
        delta = [pos2[i] - pos1[i] for i in range(3)]

        # Read strain value
        with open(join(rundir_i, f"strain-{i}"), "r") as f:
            strain_value = float(f.readline())

        optimized_geometries.append(delta)
        max_strains.append(strain_value)

    return max_strains, optimized_geometries


def plot_optimized_geometry(
    max_strains: List[float],
    optimized_geometries: List[List[float]],
    ymin: float,
    ymax: float,
    isCartesian: bool,
):
    """Plot the optimized geometry.

    :param max_strains: List of strain values.
    :param optimized_geometries: List of relative coordinates.
    :param ymin: Minimum value for the y-axis.
    :param ymax: Maximum value for the y-axis.
    :param isCartesian: coordinate type, either "lattice" - False or "cartesian" - True.
    """

    plt.rcParams.update(
        {
            "figure.figsize": (10, 7.5),
            "axes.linewidth": 4.0,
            "lines.linewidth": 3.5,
            "lines.markersize": 12.0,
            "grid.linewidth": 1.0,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "axes.edgecolor": "black",
            "axes.labelsize": 20,
            "axes.labelcolor": "black",
            "axes.axisbelow": "True",
            "legend.fontsize": 25,
            "ytick.minor.size": 6,
            "xtick.major.pad": 4,
            "ytick.major.pad": 10,
            "axes.titlesize": 30,
            "axes.titlepad": 20,
        }
    )

    fig, ax = plt.subplots()

    acoord = "cartesian" if isCartesian else "lattice"

    plt.xlabel(r"Lagrangian strain")
    plt.ylabel(r"Relative coordinate (" + acoord + ")")

    ax.xaxis.set_major_locator(
        ptk.MaxNLocator(
            nbins=6,
        )
    )
    ax.yaxis.set_major_formatter(ptk.ScalarFormatter(useOffset=True, useMathText=True))

    optimized_geometries = np.array(optimized_geometries).T

    y1, y2, y3 = optimized_geometries[:3].tolist()

    # Find reference geometry where strain is zero
    ref_index = max_strains.index(0.0)
    ref_y1, ref_y2, ref_y3 = y1[ref_index], y2[ref_index], y3[ref_index]

    plt.plot(max_strains, y1, "ro--", label=r"$\Delta 1$", zorder=3)
    plt.plot(max_strains, y2, "bs--", label=r"$\Delta 2$", zorder=2)
    plt.plot(max_strains, y3, "gd--", label=r"$\Delta 3$", zorder=1)

    plt.plot(
        [min(max_strains), max(max_strains)],
        [ref_y1, ref_y1],
        "r-",
        label=r"$\Delta$1$_{ref}$",
        zorder=3,
    )
    plt.plot(
        [min(max_strains), max(max_strains)],
        [ref_y2, ref_y2],
        "b-",
        label=r"$\Delta$2$_{ref}$",
        zorder=2,
    )
    plt.plot(
        [min(max_strains), max(max_strains)],
        [ref_y3, ref_y3],
        "g-",
        label=r"$\Delta$3$_{ref}$",
        zorder=1,
    )

    if ymin is not None and ymax is not None:
        ax.set_ylim(ymin, ymax)
    else:
        ymin = min(y1 + y2 + y3)
        ymax = max(y1 + y2 + y3)
        if abs(ymax - ymin) < 0.000000001:
            ymin -= 0.1
            ymax += 0.1
        ax.set_ylim(ymin, ymax)

    ax.grid(True, linestyle="--")

    # plt.legend(bbox_to_anchor=(1.03, 1), loc=2, borderaxespad=0., numpoints=1)
    plt.legend(loc="best")

    plt.savefig("PLOT.png", format="png", dpi=300)

    plt.close()


def main() -> None:
    parser = ArgumentParser(
        description="Python visualization tool for relaxed coordinates of atoms in the unit cell."
    )

    parser.add_argument(
        "--ATOM1",
        "--a1",
        type=int,
        help="First atom number",
        default=1,
        dest="atom1",
    )

    parser.add_argument(
        "--ATOM2",
        "--a2",
        type=int,
        help="Second atom number",
        default=2,
        dest="atom2",
    )

    parser.add_argument(
        "--y_min", type=float, dest="ymin", help="minimum value of y axis", default=None
    )

    parser.add_argument(
        "--y_max", type=float, dest="ymax", help="maximum value of y axis", default=None
    )

    args = parser.parse_args()
    atom1 = args.atom1
    atom2 = args.atom2
    ymin = args.ymin
    ymax = args.ymax

    rundir_1 = join(os.getcwd(), "rundir-1")
    if not os.path.exists(rundir_1):
        raise ValueError("Directory doesn't have rundir")

    # Get the number of atoms
    num_atoms = get_num_atoms(rundir_1)

    # Check if the specified atom numbers exceed the total number of atoms
    if atom1 > num_atoms or atom2 > num_atoms:
        raise ValueError(
            f"The specified atoms {atom1} or {atom2} exceed the total number of atoms {num_atoms}."
        )

    # Get the coordinate type
    isCartesian = is_coordinate_cartesian(rundir_1)

    # Get the optimized positions
    max_strains, optimized_geometries = get_optimized_relative_coordinates(
        atom1, atom2
    )

    # Plot the data
    plot_optimized_geometry(max_strains, optimized_geometries, ymin, ymax, isCartesian)


if __name__ == "__main__":
    main()
