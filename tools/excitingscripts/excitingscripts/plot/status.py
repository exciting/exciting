"""Python visualization tool for the RMS deviations of the SCF potential as a function of the iteration
number during the SCF loop.
"""

import os
from argparse import ArgumentParser
from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ptk
import numpy as np
from excitingtools.exciting_dict_parsers.groundstate_parser import parse_info_xml

if matplotlib.__version__.split(".")[0] == "2":
    matplotlib.style.use("classic")


def plot_status(run_dir: Union[str, Path]) -> None:
    """Python visualization tool for the RMS deviations of the SCF potential as a function
     of the iteration number during the SCF loop.

    :param run_dir: directory where exciting runs
    """
    run_dir = Path(run_dir)
    assert run_dir.exists(), "run dir not found"
    rmsd_file = run_dir / "RMSDVEFF.OUT"
    info_file = run_dir / "info.xml"

    if rmsd_file.exists():
        data = np.genfromtxt(rmsd_file)
    elif info_file.exists():
        info_xml_scl = parse_info_xml(info_file.as_posix())["scl"]
        data = [info_xml_scl[x]["rms"] for x in info_xml_scl]
    else:
        raise FileNotFoundError(f"Neither {rmsd_file} nor {info_file} found in given run_dir: {run_dir}!")

    if len(data) <= 1:
        print("Data not (yet) available for visualization.")
        return

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
        "axes.formatter.limits": (-5, 6),
    }
    plt.rcParams.update(params)

    plt.subplots_adjust(left=0.20, right=0.93, bottom=0.18, top=0.88, wspace=None, hspace=None)
    yfmt = ptk.ScalarFormatter(useOffset=True, useMathText=True)
    fig = matplotlib.pyplot.figure(figsize=(8, 5.5))
    ax = fig.add_subplot(111)

    x_label = "Iteration number"
    y_label = "RMSD $V_{SCF}$"

    ax.text(
        0.5,
        -0.17,
        x_label,
        size=font_label,
        transform=ax.transAxes,
        ha="center",
        va="center",
        rotation=0,
    )
    ax.text(
        -0.19,
        0.5,
        y_label,
        size=font_label,
        transform=ax.transAxes,
        ha="center",
        va="center",
        rotation=90,
    )

    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(2)

    plt.xticks(size=font_tick)
    plt.yticks(size=font_tick)
    plt.grid(True)

    x = np.arange(1, len(data) + 1)

    plt.plot(x, data, "b-")
    plt.plot(x, data, "go", label="calculated")

    plot_frame_fraction = 1 / 20
    x_min = x[0] - len(x) * plot_frame_fraction
    x_max = x[-1] + len(x) * plot_frame_fraction

    ax.yaxis.set_major_formatter(yfmt)
    ax.set_yscale("log")
    ax.set_xlim(x_min, x_max)
    ax.set_axisbelow(True)
    plt.legend(loc="upper right")

    plt.savefig("PLOT_status.png", orientation="portrait", format="png", dpi=300)


def main() -> None:
    parser = ArgumentParser(
        description="Python visualization tool for the RMS deviations of the SCF potential"
        " as a function of the iteration number during the SCF loop."
    )

    parser.add_argument("--run-directory", "-r",
                        default=os.getcwd(),
                        nargs=1,
                        dest="run_directory",
                        help="root path for files that are created by this script")

    parser.add_argument("-sh", "--show",
                        action="store_true",
                        help="show plot")

    args = parser.parse_args()

    plot_status(args.run_directory)

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
