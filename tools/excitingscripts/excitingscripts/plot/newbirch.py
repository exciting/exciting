"""Fit energy-vs-volume curves using the Birch-Murnaghan equation of state (**BM-EoS**) in polynomial form.

Located at `excitingscripts/plot/newbirch.py`.

Call as:

```bash
python3 -m excitingscripts.plot.newbirch
```
"""

import os
import pathlib
from argparse import ArgumentParser
from typing import Union, Tuple

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ptk
import numpy as np
import pylab as pyl
from excitingscripts.utils.utils import sort_lists_by_first_list

if matplotlib.__version__.split(".")[0] == "2":
    matplotlib.style.use("classic")


def sortstrain(s: np.ndarray, e: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    inds = s.argsort()
    s = s[inds]
    e = e[inds]
    return s, e


def main() -> None:
    parser = ArgumentParser(
        description="""Fit energy-vs-volume curves using the Birch-Murnaghan equation of state
                                        (BM-EoS) in polynomial form."""
    )

    parser.add_argument(
        "--root-directory",
        "-r",
        type=Union[str, pathlib.Path],
        default=[os.getcwd()],
        nargs=1,
        dest="root_directory",
        help="path for root directory",
    )

    parser.add_argument("-sh", "--show", action="store_true", help="show plot")

    args = parser.parse_args()

    try:
        energy_volume_data = np.loadtxt(f"{args.root_directory[0]}/energy-vs-volume")
    except FileNotFoundError:
        print("File energy-vs-volume not found!")

    strain = energy_volume_data[:, 0] ** (-2 / 3)
    energy = energy_volume_data[:, 1]
    strain, energy = sort_lists_by_first_list(strain, energy)

    bohr_radius = 0.529177
    joule2hartree = 4.3597482
    unitconv = joule2hartree / bohr_radius**3 * 10**3

    order_of_fit = 3
    fitr = np.polyfit(strain, energy, order_of_fit)
    curv = np.poly1d(fitr)
    bulk = np.poly1d(np.polyder(fitr, 2))
    bpri = np.poly1d(np.polyder(fitr, 3))
    vmin = np.roots(np.polyder(fitr))

    dmin = []
    for i in range(len(vmin)):
        if abs(vmin[i].imag) < 1.0e-10:
            if strain[0] <= vmin[i] <= strain[-1]:
                if bulk(vmin[i]) > 0:
                    dmin.append(vmin[i].real)

    chi = 0
    for i in range(len(energy)):
        chi = chi + (energy[i] - curv(strain[i])) ** 2
    chi = np.sqrt(chi) / len(energy)

    xvol = np.linspace(strain[0], strain[-1], 100)

    rvol = []
    rene = []
    for i in range(len(xvol)):
        rvol.append(xvol[i] ** (-3 / 2))
        rene.append(curv(xvol[i]))

    rstr = []
    for i in range(len(strain)):
        rstr.append(strain[i] ** (-3 / 2))

    rmin = []
    emin = []
    for i in range(len(dmin)):
        rmin.append(dmin[i] ** (-3 / 2))
        emin.append(curv(dmin[i]))

    # Plot settings

    xlabel = "Volume [Bohr\u00b3]"
    ylabel = r"Energy [Ha]"
    if os.path.exists("quantum-espresso"):
        ylabel = r"Energy [Ry]"
    if os.path.exists("vasp"):
        ylabel = r"Energy [Ry]"

    fontlabel = 20
    fonttick = 16

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
    plt.subplots_adjust(
        left=0.21, right=0.93, bottom=0.18, top=0.88, wspace=None, hspace=None
    )

    yfmt = ptk.ScalarFormatter(useOffset=True, useMathText=True)

    figure = plt.figure(figsize=(8, 5.5))
    ax = figure.add_subplot(111)
    ax.set_xlabel(xlabel, labelpad=10, fontsize=fontlabel)
    ax.set_ylabel(ylabel, labelpad=10, fontsize=fontlabel)
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(2)

    plt.xticks(size=fonttick)
    plt.yticks(size=fonttick)
    pyl.grid(True)
    plt.plot(rvol, curv(xvol), "b-", label="birch-murnaghan fit")
    plt.plot(rstr, energy, "go", label="calculated")
    plt.plot(rmin, emin, "ro")
    plt.legend(loc=9, borderaxespad=0.8, numpoints=1)

    ymax = max(max(curv(xvol)), max(energy))
    ymin = min(min(curv(xvol)), min(energy))
    dxx = abs(max(rvol) - min(rvol)) / 18
    dyy = abs(ymax - ymin) / 18
    ax.yaxis.set_major_formatter(yfmt)
    ax.set_xlim(min(rvol) - dxx, max(rvol) + dxx)
    ax.set_ylim(ymin - dyy, ymax + dyy)

    ax.xaxis.set_major_locator(ptk.MaxNLocator(7))

    plt.tight_layout()
    plt.savefig("PLOT.png", orientation="portrait", format="png", dpi=300)

    if args.show:
        plt.show()

    # Print fit results

    if len(dmin) > 1:
        print("##############################################\n")
        print("WARNING: Multiple minima are found!\n")
        print("##############################################\n")

    fmt = "%11.5f"
    amt = "%10.4f"
    bmt = "%9.3f"
    lmt = "%10.2f"

    for i in range(len(dmin)):
        x0 = dmin[len(dmin) - 1 - i]
        v0 = rmin[len(rmin) - 1 - i]
        a0sc = (1 * v0) ** (0.33333333333)
        abcc = (2 * v0) ** (0.33333333333)
        afcc = (4 * v0) ** (0.33333333333)

        derivV2 = 4 / 9 * x0**5 * bulk(x0)
        derivV3 = -20 / 9 * x0 ** (13 / 2) * bulk(x0) - 8 / 27 * x0 ** (15 / 2) * bpri(
            x0
        )
        b0 = derivV2 / x0 ** (3 / 2) * unitconv
        bp = -1 - x0 ** (-3 / 2) * derivV3 / derivV2

        print(
            "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        )
        print(
            "     V0        B0         Bp        a-sc       a-bcc      a-fcc     log(chi)"
        )
        print(
            fmt % (v0),
            bmt % (b0),
            bmt % (bp),
            amt % (a0sc),
            amt % (abcc),
            amt % (afcc),
            lmt % (np.log10(chi)),
        )
        print(
            "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        )

    if len(dmin) == 0:
        print("##############################################\n")
        print("WARNING: No minimum in the given xrange!\n")
        print("##############################################\n")


if __name__ == "__main__":
    main()
