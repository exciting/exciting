"""Python visualization tool for the cartesian components of the position of the center of mass during relaxation"""

from os.path import join
from excitingscripts.utils.utils import get_structure_optimizations_properties, get_num_atoms, is_coordinate_cartesian
from argparse import ArgumentParser
from typing import List, Tuple
import matplotlib.pyplot as plt
import matplotlib.ticker as ptk
import numpy as np


def plot_center_of_mass(center_of_mass_data: List[Tuple[float, float, float]], run_dir: str, show: bool, dpi: int):
    """ Plot the center of mass components, then save the plot to a file.

    :param center_of_mass_data: A list of tuples containing the center of mass data.
    :param run_dir: directory where exciting runs
    :param show: Whether to display the plot.
    :param dpi: Resolution in DPI for the saved plot.
    """
    plt.rcParams.update({
        'figure.figsize': (10, 7.5),
        'axes.linewidth': 4.0,
        'lines.markersize': 10,
        'lines.linewidth': 2.0,
        'grid.linewidth': 1.0,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'ytick.minor.size': 6,
        'xtick.major.pad': 10,
        'ytick.major.pad': 10,
        'axes.formatter.limits': (-5, 6),
        'axes.labelsize': 20,
        'axes.axisbelow': 'True',
        'legend.fontsize': 25,
    })

    fig, ax = plt.subplots()

    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(2)

    ax.set_xlabel(u'Optimization steps')
    ax.set_ylabel(r'Center of mass (cartesian)')
    ax.grid(True, linestyle='--')

    x = list(range(len(center_of_mass_data)))
    r_x = [abs(t[0]) for t in center_of_mass_data]
    r_y = [abs(t[1]) for t in center_of_mass_data]
    r_z = [abs(t[2]) for t in center_of_mass_data]

    ax.plot(x, r_x, 'ro-', label=u'R$_x$')
    ax.plot(x, r_y, 'bo-', label=u'R$_y$')
    ax.plot(x, r_z, 'go-', label=u'R$_z$')

    ax.xaxis.set_major_locator(ptk.MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(ptk.ScalarFormatter(useOffset=True, useMathText=True))

    # Calculate x-axis limits
    iter_count = len(center_of_mass_data)
    xmin = 0 - (iter_count - 1) / 20.0
    xmax = (iter_count - 1) + (iter_count - 1) / 20.0

    if iter_count == 1:
        xmin = -1
        xmax = 1

    ax.set_xlim(xmin, xmax)

    # Calculate y-axis limits
    ymin = min(min(r_x), min(r_y), min(r_z))
    ymax = max(max(r_x), max(r_y), max(r_z))

    if abs(ymax - ymin) < 1.e-9:
        ymin -= 1
        ymax += 1
    else:
        dy = abs(ymax - ymin) / 18.0
        ymin = ymin - dy
        ymax = ymax + dy

    ax.set_ylim(ymin, ymax)

    # Adjust legend
    plt.legend(borderaxespad=.8)

    plt.savefig(join(run_dir, "PLOT.png"), format='png', dpi=dpi)

    if show:
        plt.show()

    plt.close()


def main() -> None:
    parser = ArgumentParser(
        description="Python visualization tool for the cartesian components of the position of the center of mass during relaxation")

    # Define command line arguments
    parser.add_argument("directory",
                        type=str,
                        help="Directory where exciting runs")

    parser.add_argument("-sh", "--show",
                        action="store_true",
                        help="Show plot in a window")

    parser.add_argument("--dpi",
                        type=int,
                        default=300,
                        help="Resolution in DPI for saved plot")

    args = parser.parse_args()
    run_dir = args.directory

    # Extract center of mass data from the file
    try:
        center_of_mass_data = get_structure_optimizations_properties(run_dir, "Center of mass")
    except ValueError:
        # Probably the COM data doesn't exist, but we can work with Atomic positions
        atomic_positions = get_structure_optimizations_properties(run_dir, "Atomic positions")
        center_of_mass_data = []
        for atomic_position in atomic_positions:
            pos_matrix = []
            for i in list(atomic_position.keys()):
                pos_matrix.append(atomic_position[i])
            com = np.mean(pos_matrix, axis=0).tolist()
            center_of_mass_data.append(com)

    # Plot center of mass components
    plot_center_of_mass(center_of_mass_data, run_dir, args.show, args.dpi)


if __name__ == "__main__":
    main()

