import os
import pathlib
import matplotlib as mpl
import matplotlib.pyplot as plt

def get_exciting_root() -> str:
    """Get exciting's root directory.

    :return exciting root directory.
    """
    if not 'EXCITINGROOT' in os.environ.keys():
        print(f"wd in get_exciting_root: {os.getcwd()}")
        excitingroot = str(pathlib.Path(os.getcwd()).parents[3]) 
        os.environ['EXCITINGROOT'] = excitingroot
    else:
        excitingroot = os.environ['EXCITINGROOT']
    assert "exciting_smp" in os.listdir(os.path.join(excitingroot, "bin")), f"Exciting binary 'exciting_smp' was not found at exciting root path: {excitingroot}"
    return excitingroot

# Modified from C. Vorwerk,
# https://git.physik.hu-berlin.de/sol/toolbox/-/blob/master/notebooks/BSE/Example.ipynb
def set_plot_parameters(figure_width: float, figure_height: float,
                        fontsize: float):
    """Sets 'seaborn-notebook' as plt.style and defines plot parameters
    such as fontsize and grids. 
    
    :param figure_width (float): Figure width in inches
    :param figure_height (float): Figure height in inches
    :param fontsize (float): Font size in pt


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
    plt.rcParams['legend.fontsize'] = 0.75*plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = 0.75*plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = 0.75*plt.rcParams['font.size']
    plt.rcParams['lines.linewidth']=3


