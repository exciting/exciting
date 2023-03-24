import os
import pathlib
from typing import Union
import matplotlib as mpl
import matplotlib.pyplot as plt
import re
from nbformat import read as read_nb


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


def re_input(nb_path: Union[str, pathlib.Path], title: str) -> str:
    """ Extract an exciting input file string from a Jupyter notebook file.

    Only return the first match. The developer should therefore
    choose a unique identifying title, for example:

      <title>DiamondForTutorial</title>

    for the input they wish to parse.

    :param nb_path: Notebook file location.
    :param title: Title string between the XML tags.
    :return parsed input.xml.
    """
    # This can likely be improved
    pattern = f'(<input>[^;]*<title>{title}[^;]*<[/]input>)'
    nb = read_nb(nb_path, as_version=4)
    for cell in filter(lambda x: x["cell_type"] == "markdown", nb["cells"]):
        matches = re.findall(pattern, cell["source"], flags = re.MULTILINE)
        if len(matches) != 0:
            return matches[0]
    return 'NULL'

def get_input_xml_from_notebook(nb_path: Union[str, pathlib.Path], ident: str) -> str:
    """ Extract an exciting input file string from a Jupyter notebook file.

    Only return the first match. The input is identified by a (hidden) html tag, 
    where the html class name identifies the cell to parse.

    The structure of the cells must be:

        <span class="{ident}"></span>
        ```xml
        <content/>
        ```

    where `{ident}` is replaced with the same string as in this function.

    :param nb_path: Notebook file location (excluding file suffix).
    :param ident: Identifier of input string cell.
    :return parsed input.xml.
    """
    nb = read_nb(f"{nb_path}.ipynb", as_version=4)
    # check markdown cells
    for cell in filter(lambda x: x["cell_type"] == "markdown", nb["cells"]):
        if ident in cell["source"]:
            cell_content = cell["source"]
            break
    # delete first two (containing ident and markdown formatting) and last (contain formatting) line
    cell_content = "\n".join(cell_content.split("\n")[2:-1])
    return cell_content