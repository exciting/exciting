"""This script allows for the visualization of the force-vs-displacement curve."""

import json
from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ptk

from excitingscripts.utils.utils import sort_lists_by_first_list

if matplotlib.__version__.split(".")[0] == "2":
    matplotlib.style.use("classic")


def plot_atomforce() -> None:
    """This script allows for the visualization of the force-vs-displacement curve."""
    phonon_results_file = Path("phonon_results.json")
    if not phonon_results_file.exists():
        raise FileNotFoundError("file phonon_results.json not found")

    forces = []
    displacements = []

    with open(phonon_results_file) as fid:
        phonon_data = json.load(fid)["results"]

    for displacement_string, result in phonon_data.items():
        displacements.append(float(displacement_string))
        forces.append(result["force"])

    displacements, forces = sort_lists_by_first_list(displacements, forces)

    # calculating plot ranges
    plot_frame_fraction = 1 / 18
    max_force = max(forces)
    min_force = min(forces)
    x_increment = abs(displacements[-1] - displacements[0]) * plot_frame_fraction
    y_increment = abs(max_force - min_force) * plot_frame_fraction

    x_min = displacements[0] - x_increment
    x_max = displacements[-1] + x_increment
    y_min = min_force - y_increment
    y_max = max_force + y_increment

    # set defaults parameters for the plot
    font_label = 18
    font_tick = 14

    params = {
        "ytick.minor.size": 6,
        "xtick.major.pad": 8,
        "ytick.major.pad": 4,
        "patch.linewidth": 2.0,
        "axes.linewidth": 2.0,
        "lines.linewidth": 1.8,
        "lines.markersize": 8.0,
        "axes.formatter.limits": (-4, 6),
    }

    plt.rcParams.update(params)
    plt.subplots_adjust(left=0.21, right=0.93, bottom=0.18, top=0.88, wspace=None, hspace=None)

    yfmt = ptk.ScalarFormatter(useOffset=True, useMathText=True)
    fig = plt.figure(figsize=(8, 5.5))
    ax = fig.add_subplot(111)
    y_label = "Atomic force [Ha/Bohr]"
    x_label = "Displacement $u$ [alat]"

    ax.text(0.5, -0.17, x_label, size=font_label, transform=ax.transAxes, ha="center", va="center", rotation=0)
    ax.text(-0.19, 0.5, y_label, size=font_label, transform=ax.transAxes, ha="center", va="center", rotation=90)

    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(2)

    plt.xticks(size=font_tick)
    plt.yticks(size=font_tick)
    plt.grid(True)

    plt.plot(displacements, forces, "r-")
    plt.plot(displacements, forces, "go", label="calculated")

    plt.legend(borderaxespad=0.8, numpoints=1)

    ax.yaxis.set_major_formatter(yfmt)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.xaxis.set_major_locator(ptk.MaxNLocator(6))
    ax.yaxis.set_major_locator(ptk.MaxNLocator(5))
    ax.set_axisbelow(True)

    plt.savefig("PLOT_atomforce.png", orientation="portrait", format="png", dpi=300)

def main() -> None:
    """Parser and function call."""
    parser = ArgumentParser(description="Plot force-vs-displacement curves.")

    parser.add_argument("-sh", "--show",
                        action="store_true",
                        help="show plot")

    args = parser.parse_args()

    plot_atomforce()

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
