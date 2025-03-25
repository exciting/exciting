"""Python visualization tool for the total torque during relaxation."""

from os.path import join
import numpy as np
from argparse import ArgumentParser
from typing import List, Dict, Tuple
from excitingscripts.utils.utils import get_structure_optimizations_properties
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ptk

# Use classic style if Matplotlib version is 2
if matplotlib.__version__.split(".")[0] == "2":
    matplotlib.style.use('classic')


def plot_torque(torque_data: List[Tuple[float, float, float]], run_dir: str, show: bool, dpi: int, tol=1.0e-6):
    """ Plot the torque components and magnitude, then save the plot to a file.

    :param torque_data: list of torque data to plot.
    :param run_dir: directory where exciting runs.
    :param show: whether to display the plot.
    :param dpi: resolution in DPI for the saved plot.
    :param tol: determines the lowest value possible
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

    ax.set_xlabel('Optimization steps')
    ax.set_ylabel('|Total torque| [Ha]')
    ax.grid(True, linestyle='--')

    torque_data = np.absolute(torque_data)

    x = list(range(len(torque_data)))
    u_x = [t[0] for t in torque_data]
    u_y = [t[1] for t in torque_data]
    u_z = [t[2] for t in torque_data]

    for i in range(len(u_x)):
        if u_x[i] < tol:
            u_x[i] = tol
        if u_y[i] < tol:
            u_y[i] = tol
        if u_z[i] < tol:
            u_z[i] = tol

    ax.plot(x, u_x, 'ro-', label=u'T$_x$')
    ax.plot(x, u_y, 'bo-', label=u'T$_y$')
    ax.plot(x, u_z, 'go-', label=u'T$_z$')

    ax.xaxis.set_major_locator(ptk.MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(ptk.ScalarFormatter(useOffset=True, useMathText=True))

    # Calculate x-axis limits
    iter_count = len(torque_data)
    xmin = 0 - (iter_count - 1) / 20.0
    xmax = (iter_count - 1) + (iter_count - 1) / 20.0

    if iter_count == 1:
        xmin = -1
        xmax = 1

    ax.set_xlim(xmin, xmax)

    # Calculate y-axis limits
    ymin = min(min(u_x), min(u_y), min(u_z))
    ymax = max(max(u_x), max(u_y), max(u_z))

    if abs(ymax - ymin) < 1.e-9:
        ymin -= 1
        ymax += 1
        ax.set_ylim(ymin, ymax)
    else:
        ax.set_yscale('log')

    # Adjust legend
    plt.legend(borderaxespad=.8)

    plt.savefig(join(run_dir, "PLOT.png"), format='png', dpi=dpi)

    if show:
        plt.show()

    plt.close()


def main() -> None:
    parser = ArgumentParser(
        description="Python visualization tool for the total torque during relaxation.")

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

    # Extract torque data from the file
    torque_data = get_structure_optimizations_properties(run_dir, 'Total torque')

    # Plot torque components and magnitudes
    plot_torque(torque_data, run_dir, args.show, args.dpi)


if __name__ == "__main__":
    main()

