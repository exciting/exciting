import pathlib
import re
from argparse import ArgumentParser
from typing import Union

import matplotlib.pyplot as plt
import numpy as np

# Plot settings
figcolor = 'white'
dpi = 300
fig = plt.figure(figsize=(14.5,10),dpi=dpi)
fig.patch.set_edgecolor(figcolor)
fig.patch.set_facecolor(figcolor)

plt.rcParams['axes.linewidth' ] = 4.0
plt.rcParams['xtick.labelsize'] = 30
plt.rcParams['ytick.labelsize'] = 30
plt.rcParams['axes.edgecolor' ] = 'black'
plt.rcParams['axes.labelsize' ] = 40
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['axes.axisbelow' ] = 'True'
plt.rcParams['legend.fontsize'] = 40
plt.rcParams['xtick.major.pad'] = 10
plt.rcParams['ytick.major.pad'] = 10

colors=['firebrick','mediumblue','g','y','c','m','k']

ax1 = fig.add_subplot(111)

ax1.xaxis.set_label_position('bottom')
ax1.set_xlabel('Energy [eV]', labelpad=19)
ax1.set_ylabel(r'Im $\epsilon_M$', labelpad=13)

ax1.axhline(y=0, linestyle="dashed", linewidth=3, color="black")

for line in ax1.get_xticklines() + ax1.get_yticklines():
    line.set_markersize(10)
    line.set_markeredgewidth(2)

def plot_spectra(plot_file_path : Union[str, pathlib.Path], color_index) -> None:
    """Plot imaginary part of the macroscopic dielectric function.

    :param plot_file_path: Path to file containing data wanted for plot.
    :param color_index: Index needing for plots with different colors for each calculation.
    """
    file_name = plot_file_path.split("/")[-1]
    file_name = re.split('_|-', file_name)

    legend = None
    if "BSE" in file_name: legend = "BSE "
    if "FXCRPA" in file_name: legend = "RPA "
    if "FXCALDA" in file_name: legend = "ALDA "
    if "FXCLRCstatic" in file_name: legend = "LRCstatic "
    if "FXCRBO" in file_name: legend = "RBO "
    if "FXCLRCdyn" in file_name: legend = "LRCdyn "
    if "FXCMB1" in file_name: legend = "MB1 "
    if "TDA" in file_name: legend = legend + "TDA "
    if  not "TDA" in file_name and "BSE" in file_name: legend = legend + "full "
    if "singlet" in file_name: legend = legend + "singlet"
    if "triplet" in file_name: legend = legend + "triplet"
    if "RPA" in file_name: legend = legend + "RPA"
    if "IP" in file_name: legend = legend + "IP"
    if  not "NLF" in file_name and not "BSE" in file_name: legend = legend + "(LFE) "
    if "NLF" in file_name: legend = legend + "(no-LFE) "
    if file_name[-2][0:2]=="OC" and "LOSS" in file_name: legend = legend + "Optical(%s)"%(file_name[-2][2:])

    spectrum_data = np.genfromtxt(plot_file_path)
    ax1.plot(spectrum_data[:, 0], spectrum_data[:, 2], color=colors[color_index], label=legend, lw=3.0)

def main() -> None:
    parser = ArgumentParser(description="Plot imaginary part of the macroscopic dielectric function.")

    parser.add_argument("--plot-files", "-f",
                        nargs='+',
                        dest="plotfiles",
                        help="names of the files containg plot data")

    parser.add_argument("-sh", "--show",
                        action="store_true",
                        help="show plot")

    parser.add_argument("-pn", "--plot_name",
                        type=str,
                        default="PLOT",
                        help="plot name")

    args = parser.parse_args()

    for index, plot_file in enumerate(args.plotfiles):
        plot_spectra(plot_file, index)

    for line in ax1.get_xticklines() + ax1.get_yticklines():
        line.set_markersize(8)
        line.set_markeredgewidth(3)

    leg = ax1.legend(loc="upper left", borderaxespad=0.5, numpoints=1, fontsize=30)
    leg.get_frame().set_linewidth(4.0)
    leg.get_frame().set_edgecolor("grey")
    leg.draw_frame(True)

    plt.tight_layout()
    plt.savefig(f"{args.plot_name}.png", orientation='portrait', format='png', dpi=300)

    if args.show:
        plt.show()

if __name__ == "__main__":
    main()