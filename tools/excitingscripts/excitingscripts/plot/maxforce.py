"""Python visualization tool for the maximum amplitude of the force on the atoms during relaxation."""

import math
import os
from argparse import ArgumentParser
from os.path import join
from typing import List

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ptk
import numpy as np
from excitingscripts.utils.utils import get_structure_optimizations_properties

if matplotlib.__version__.split(".")[0] == "2":
    matplotlib.style.use('classic')


def plot_forces(forces: List[List], run_dir: str, show: bool, dpi: int):
    """ Plot the torque components and magnitude, then save the plot to a file.

    :param forces: list of forces and target.
    :param run_dir: directory where exciting runs.
    :param show: whether to display the plot.
    :param dpi: resolution in DPI for the saved plot.
    """
    # Plot settings
    plt.rcParams.update({
        'figure.figsize': (10, 7.5),
        'axes.linewidth': 4.0,
        'lines.markersize': 10,
        'lines.linewidth': 3.0,
        'grid.linewidth': 1.0,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'axes.edgecolor': 'black',
        'axes.labelsize': 20,
        'axes.labelcolor': 'black',
        'axes.axisbelow': 'True',
        'legend.fontsize': 25,
        'ytick.minor.size': 6,
        'xtick.major.pad': 10,
        'ytick.major.pad': 10,
        'axes.titlesize': 30,
        'axes.titlepad': 20
    })
    # Create figure and axis
    fig, ax = plt.subplots()

    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(2)

    plt.xticks()
    plt.yticks()

    forces = np.array(forces)
    max_forces = forces[:, 0].tolist()
    goal = forces[0, 1]

    # Plot the maximum forces and the goal
    ax.plot([0, len(max_forces) - 1], [goal, goal], 'b-', label='Goal', )
    ax.plot(max_forces, 'ro-', label='Run', )

    # Set labels
    ax.set_xlabel('Optimization steps')
    ax.set_ylabel('Maximum force [Ha/Bohr]')
    ax.grid(True, linestyle='--')

    # Ensure that only integer ticks are shown on the x-axis
    ax.xaxis.set_major_locator(ptk.MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(ptk.ScalarFormatter(useOffset=True, useMathText=True))

    # Setting the y-axis limits between order(goal for convergence) - 1 order(initial max_force) + 1
    order_upper = math.floor(max(math.log(goal, 10), math.log(max_forces[0], 10)))
    order_bottom = math.ceil(min(math.log(goal, 10), math.log(max_forces[-1], 10)))
    ax.set_ylim(10 ** (order_bottom - 1), 10 ** (order_upper + 1.5))
    ax.set_yscale('log')

    # Add legend
    ax.legend(borderaxespad=0.8)

    # Save the plot as a PNG file
    plt.savefig(join(run_dir, "PLOT.png"), format='png', dpi=dpi)

    # Show the plot if the show flag is set
    if show:
        plt.show()

    plt.close()


def main() -> None:
    parser = ArgumentParser(
        description="Python visualization tool for the maximum amplitude of the force on the atoms during relaxation.")

    parser.add_argument("--run-directory", "-r",
                        default=os.getcwd(),
                        nargs=1,
                        dest="run_directory",
                        help="root path for files that are created by this script")

    parser.add_argument("-sh", "--show",
                        action="store_true",
                        help="Show plot in a window")

    parser.add_argument("--dpi",
                        type=int,
                        default=300,
                        help="Resolution in DPI for saved plot")

    args = parser.parse_args()
    run_dir = args.run_directory[0]

    forces = get_structure_optimizations_properties(run_dir, "Maximum force")

    plot_forces(forces, run_dir, args.show, args.dpi)


if __name__ == "__main__":
    main()

