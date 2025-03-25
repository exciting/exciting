"""Plot derivatives.

This is a very important tool that allows to represent the dependence of the calculated derivatives of the
energy-vs-displacement and force-vs-displacement curves on

 * the range of points included in the fitting procedure ("maximum displacement u"),
 * the maximum degree of the polynomial used in the fitting procedure ("n").
"""

import json
from argparse import ArgumentParser
from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ptk
from excitingscripts.utils.utils import (
    get_decimal_decomposition,
    get_prettified_scientific_notation,
)

if matplotlib.__version__.split(".")[0] == "2":
    matplotlib.style.use("classic")


def find_first_none(inp: list) -> Union[int, None]:
    """Finds the first none in the list.

    :param inp: input list
    :return: the index of the first none in the inp list
    """
    for i, entry in enumerate(inp):
        if entry is None:
            return i
    return None


def plot_checkderiv(
    quantity: str = "energy",
    y_min_arg: Union[float, None] = None,
    y_max_arg: Union[float, None] = None,
) -> None:
    """Plot the derivative of the energy-vs-displacement or force-vs-displacement curves.

    :param quantity: of interest, could be 'energy' or 'force' or 'strain'
    :param y_min_arg: lower limit of the y-axis for plotting
    :param y_max_arg: upper limit of the y-axis for plotting
    """
    if Path("energy-vs-strain").exists():
        unit = "GPa"
        x_label = "Maximum lagrangian strain"
        if Path("planar").exists():
            unit = "N/m"
    elif Path("phonon_results.json").exists():
        unit = r"cm$^-\!$" + "\u00b9"
        x_label = "Maximum displacement $u$ [alat]"
    else:
        raise ValueError("No tutorial file found.")

    input_file = Path(f"checkfit_{quantity}_results.json")
    if not input_file.exists():
        raise FileNotFoundError(f"File {input_file} not found.")
    with open(input_file) as fid:
        file_content = json.load(fid)

    startorder = 0
    if Path("startorder").exists():
        startorder = int(open("startorder", "r").readline().strip().split()[0])

    order = file_content["order_of_derivative"]
    ylabel = f"Derivative of order ${order}$"
    full_order = order + startorder

    y1 = []
    y2 = []
    y3 = []
    x_values = []  # could be strain or displacements

    if Path("energy-vs-strain").exists():
        for data in file_content["fits"]:
            x_values.append(data["max_strain"])
            frequencies = data["derivatives"]

            try:
                y1.append(frequencies[0]["value"])
            except:
                y1.append(None)

            try:
                y2.append(frequencies[2]["value"])
            except:
                y2.append(None)

            try:
                y3.append(frequencies[4]["value"])
            except:
                y3.append(None)
    else:
        for data in file_content["fits"]:
            x_values.append(data["max_displacement"])  # TODO: other key here?
            frequencies = data["frequencies"]
            y1.append(frequencies[0])
            y2.append(frequencies[2])
            y3.append(frequencies[4])

    y1 = y1[: find_first_none(y1)]
    y2 = y2[: find_first_none(y2)]
    y3 = y3[: find_first_none(y3)]

    # manipulate data for a better plot
    total_y_min = min(y1 + y2 + y3)
    exponent_factor = 10 ** get_decimal_decomposition(total_y_min)[1]

    y1 = [(x - total_y_min) / exponent_factor for x in y1]
    y2 = [(x - total_y_min) / exponent_factor for x in y2]
    y3 = [(x - total_y_min) / exponent_factor for x in y3]

    y_min = min(y1 + y2 + y3)
    y_max = max(y1 + y2 + y3)
    plot_frame_fraction = 1 / 18
    dyy = abs(y_max - y_min) * plot_frame_fraction
    y_min = y_min - dyy
    y_max = y_max + dyy

    if y_min_arg is not None:
        y_min = (float(y_min_arg) - total_y_min) / exponent_factor
    if y_max_arg is not None:
        y_max = (float(y_max_arg) - total_y_min) / exponent_factor

    x_min = min(x_values)
    x_max = max(x_values)
    dxx = abs(x_max - x_min) * plot_frame_fraction
    x_min -= dxx
    x_max += dxx

    # set defauls parameters for the plot
    fontlabel = 18
    fonttick = 14
    fonttext = 18
    fontlimits = 12

    params = {
        "ytick.minor.size": 6,
        "xtick.major.pad": 8,
        "ytick.major.pad": 4,
        "patch.linewidth": 2.0,
        "axes.linewidth": 2.0,
        "lines.linewidth": 3.5,
        "lines.markersize": 12.0,
        "mathtext.fontset": "stixsans",
        "axes.formatter.limits": (-8, 8),
    }

    plt.rcParams.update(params)
    plt.subplots_adjust(
        left=0.20, right=0.78, bottom=0.18, top=0.88, wspace=None, hspace=None
    )

    yfmt = ptk.ScalarFormatter(useOffset=True, useMathText=True)
    fig = matplotlib.pyplot.figure(figsize=(8, 5.5))
    ax = fig.add_subplot(111)

    ax.text(
        0.5,
        -0.17,
        x_label,
        size=fontlabel,
        transform=ax.transAxes,
        ha="center",
        va="center",
        rotation=0,
    )
    ax.text(
        0.0,
        1.05,
        get_prettified_scientific_notation(total_y_min, unit),
        size=fonttext,
        color="#00008B",
        transform=ax.transAxes,
        ha="left",
        va="center",
        rotation=0,
    )

    ax.text(
        0.62,
        1.033,
        "m = "
        + get_prettified_scientific_notation(y_min * exponent_factor + total_y_min),
        size=fontlimits,
        color="#00008B",
        transform=ax.transAxes,
        ha="left",
        va="center",
        rotation=0,
    )
    ax.text(
        0.62,
        1.080,
        "M = "
        + get_prettified_scientific_notation(y_max * exponent_factor + total_y_min),
        size=fontlimits,
        color="#00008B",
        transform=ax.transAxes,
        ha="left",
        va="center",
        rotation=0,
    )

    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(2)

    plt.xticks(size=fonttick)
    plt.yticks(size=fonttick)
    plt.grid(True)
    plt.ylabel(ylabel, size=fontlabel)

    if len(y1) > 0:
        plt.plot(x_values[: len(y1)], y1, "ro-", label=f"n={full_order}")
    if len(y2) > 0:
        plt.plot(x_values[: len(y2)], y2, "bo-", label=f"n={full_order + 2}")
    if len(y3) > 0:
        plt.plot(x_values[: len(y3)], y3, "go-", label=f"n={full_order + 4}")

    plt.legend(bbox_to_anchor=(1.03, 1), loc=2, borderaxespad=0.0)

    ax.yaxis.set_major_formatter(yfmt)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    ax.xaxis.set_major_locator(ptk.MaxNLocator(6))
    ax.set_axisbelow(True)

    plt.savefig("PLOT_checkderiv.png", orientation="portrait", format="png", dpi=300)


def main() -> None:
    """Parser and function call."""
    parser = ArgumentParser(
        description="""Plot the derivative of the energy-vs-displacement or force-vs-displacement 
                                        curves."""
    )

    parser.add_argument(
        "quantity",
        type=str,
        help="which quantity to plot, energy, force, or strain.",
        default="energy",
    )
    parser.add_argument(
        "--y_min", type=float, help="minimum value of y axis", default=None
    )
    parser.add_argument(
        "--y_max", type=float, help="maximum value of y axis", default=None
    )

    parser.add_argument("-sh", "--show", action="store_true", help="show plot")

    args = parser.parse_args()

    plot_checkderiv(args.quantity, args.y_min, args.y_max)

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
