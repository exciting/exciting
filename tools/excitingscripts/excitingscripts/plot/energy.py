"""Visualize energy-vs-strain curves.

Located at `excitingscripts/plot/energy.py`. 

Call as:

```bash
python3 -m excitingscripts.plot.energy
```
"""

import json
import os
import sys
from argparse import ArgumentParser
from typing import List, Tuple

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ptk
import numpy as np
import pylab as pyl

if matplotlib.__version__.split(".")[0] == "2":
    matplotlib.style.use('classic')

def sortstrain(strain: List[float], energy: List[float]) -> Tuple[List[float], List[float]]:
    """Sort strain values and also sort energy values based on the index of the sorted strain list.

    :param strain: List containing strain values.
    :param strain: List containing energy values.
    :returns: Lists containing sorted strain and energy values.
    """
    sorted_strain = sorted(strain)
    sorted_energy = [energy[i] for i in np.argsort(np.array(strain))]

    return sorted_strain, sorted_energy

def main() -> None:
    parser = ArgumentParser(description="Plot energy curves.")

    parser.add_argument("-sh", "--show",
                        action="store_true",
                        help="show plot")

    args = parser.parse_args()

    ylabel = r'Energy [Ha]'

    inpf = None

    if os.path.exists('energy-vs-strain'):
        inpf = 'energy-vs-strain'
        xlabel = r'Lagrangian strain'

    if os.path.exists('energy-vs-displacement'):
        inpf = 'energy-vs-displacement'
        xlabel = r'Displacement $u$ [alat]'

    if os.path.exists('energy-vs-volume'):
        inpf = 'energy-vs-volume'
        xlabel = 'Volume [Bohr\u00B3]'

    if os.path.exists('energy-vs-ngrid'):
        inpf = 'energy-vs-ngridk'
        xlabel = r'ngridk'

    if os.path.exists('energy-vs-alat'):
        inpf = 'energy-vs-alat'
        xlabel = r'Lattice parameter [Bohr]'

    if os.path.exists('energy-vs-step'):
        inpf = 'energy-vs-step'
        xlabel = r'Step'

    if os.path.exists('phonon_results.json'):
        inpf = 'phonon_results.json'
        xlabel = r'Displacement $u$ [alat]'

    if not inpf:
        sys.exit("\nERROR: file " + inpf + " not found!\n")

    x = []
    y = []

    if inpf == 'phonon_results.json':
        with open(inpf) as fid:
            phonon_data: dict = json.load(fid)["results"]

        for displacement_string, result in phonon_data.items():
            x.append(float(displacement_string))
            y.append(result["energy"])
    else:
        lines = np.genfromtxt(inpf)
        for index, line in enumerate(lines):
            x.append(line[0])
            y.append(line[1])

    x_sorted, y_sorted = sortstrain(x, y)

    if os.path.exists('rundir-oo'):
        xlabel = r'Interlayer distance [Bohr]'
        with open("./rundir-oo/TOTENERGY.OUT", "r") as infty_file:
            energy_zero = float(infty_file.readlines()[-2].strip())

        y_sorted.pop(-1)
        x_sorted.pop(-1)

        ndigits = 8
        with open("normalized-energy", "w") as norm_file:
            for index, value_y_sorted in enumerate(y_sorted):
                value_y_sorted = round(value_y_sorted - energy_zero, ndigits)
                y_sorted[index] = value_y_sorted
                norm_file_data = np.vstack((x_sorted[index], value_y_sorted)).T
                np.savetxt(norm_file, norm_file_data, fmt="%15.8f %15.10f")

    rmin = min(y_sorted)
    srmin = " "

    if not (os.path.exists('rundir-oo')):
        srmin = '\u2013 ' + str(abs(rmin))
        if rmin > 0: srmin = '+ ' + str(rmin)
        for i in range(len(y_sorted)): y_sorted[i] = (y_sorted[i] - rmin)

    x_step = 18
    y_step = 17

    dxx = abs(max(x_sorted) - min(x_sorted)) / x_step
    dyy = abs(max(y_sorted) - min(y_sorted)) / y_step

    xmin = min(x_sorted) - dxx
    xmax = max(x_sorted) + dxx
    ymin = min(y_sorted) - dyy
    ymax = max(y_sorted) + dyy

    if len(x_sorted) == 1:
        xmin = min(x_sorted) - dxx - 1
        xmax = max(x_sorted) + dxx + 1
        ymin = min(y_sorted) - dyy - 1
        ymax = max(y_sorted) + dyy + 1

    # Set default parameters for the plot

    fontlabel = 20
    fonttick = 16
    fonttext = 14

    params = {'ytick.minor.size': 6,
              'xtick.major.pad': 8,
              'ytick.major.pad': 4,
              'patch.linewidth': 2.,
              'axes.linewidth': 2.,
              'lines.linewidth': 1.8,
              'lines.markersize': 8.0,
              'axes.formatter.limits': (-4, 6)}

    plt.rcParams.update(params)
    plt.subplots_adjust(left=0.21, right=0.93,
                        bottom=0.18, top=0.88,
                        wspace=None, hspace=None)

    yfmt = ptk.ScalarFormatter(useOffset=True, useMathText=True)

    fig = matplotlib.pyplot.figure(figsize=(8, 5.5))

    ax = fig.add_subplot(111)

    if os.path.exists('rundir-oo'):
        ax.text(0.5, -0.17, xlabel, size=fontlabel,
                transform=ax.transAxes, ha='center', va='center', rotation=0)
    else:
        ax.text(0.5, -0.17, xlabel, size=fontlabel,
                transform=ax.transAxes, ha='center', va='center', rotation=0)
        ax.text(0.11, 1.03, srmin, size=fonttext,
                transform=ax.transAxes, ha='left', va='center', rotation=0)

    ax.text(-0.23, 0.5, ylabel, size=fontlabel,
            transform=ax.transAxes, ha='center', va='center', rotation=90)

    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(6)
        line.set_markeredgewidth(2)

    plt.xticks(size=fonttick)
    plt.yticks(size=fonttick)
    pyl.grid(True)

    plt.plot(x_sorted, y_sorted, 'r-')
    plt.plot(x_sorted, y_sorted, 'go', label='calculated')

    plt.legend(loc=9, borderaxespad=.8, numpoints=1)

    ax.yaxis.set_major_formatter(yfmt)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    ax.xaxis.set_major_locator(ptk.MaxNLocator(7))

    if inpf == 'phonon_results.json':
        ax.xaxis.set_major_locator(ptk.MaxNLocator(6))

    ax.set_axisbelow(True)

    plt.savefig('PLOT.png', orientation='portrait', format='png', dpi=300)

    if args.show:
        plt.show()

if __name__ == "__main__":
    main()
