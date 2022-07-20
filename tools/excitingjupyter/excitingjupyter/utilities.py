import os
import pathlib
import matplotlib as mpl
import matplotlib.pyplot as plt


def get_exciting_root() -> str:
    """Get exciting's root directory.

    :return exciting root directory.
    """
    if not 'EXCITINGROOT' in os.environ.keys():
        excitingroot = str(pathlib.Path(os.getcwd()).parents[3])
        os.environ['EXCITINGROOT'] = excitingroot
    else:
        excitingroot = os.environ['EXCITINGROOT']
    return excitingroot


def check_for_binary(excitingroot: str, binary='exciting_smp'):
    """ Check exciting binary is present in bin.

    :param excitingroot: exciting root directory.
    :param binary: Optional binary name.
    """
    bin = os.path.join(excitingroot, "bin")
    if binary not in os.listdir(bin):
        raise FileNotFoundError(f"Exciting binary '{binary}' was not found at exciting root path: {bin}")


def set_plot_parameters(figure_width: float, figure_height: float, fontsize: float):
    """Sets 'seaborn-notebook' as plt.style and defines plot parameters
    such as fontsize and grids.

    Modified from C. Vorwerk:
    https://git.physik.hu-berlin.de/sol/toolbox/-/blob/master/notebooks/BSE/Example.ipynb
    
    :param figure_width: Figure width in inches
    :param figure_height: Figure height in inches
    :param fontsize: Font size in pt
    """
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.style.use(['seaborn-notebook'])
    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.color'] = "grey"
    plt.rcParams['grid.linestyle'] = "-"
    plt.rcParams['grid.alpha'] = "0.3"
    plt.rcParams['figure.figsize'] = (figure_width, figure_height)
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = 0.75 * plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = 0.75 * plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = 0.75 * plt.rcParams['font.size']
    plt.rcParams['lines.linewidth'] = 3
