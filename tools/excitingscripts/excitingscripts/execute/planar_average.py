"""Extract planar-averaged electrostatic potential in a given direction.

Located at `excitingscripts/execute/planar_average.py`.

Call as:

```bash
python3 -m excitingscripts.execute.planar_average direction
```
Where <code>direction</code> is the direction along which the plane-averaged potential will be visualized.
"""

import pathlib
import xml.etree.ElementTree as ET
from argparse import ArgumentParser
from typing import Union

import numpy as np


def execute_planar_average(potential_file: Union[str, pathlib.Path], direction: str) -> None:
    """Extract planar-averaged electrostatic potential in a given direction.
    Script included in "execute" directory for consistency, due to classification of the old "tutorial scripts"

    :param potential_file: File containing electrostatic potential data.
    :param direction: Direction along which potential needs to be averaged.
    """
    direction = direction.upper()
    if direction not in ["X", "Y", "Z"]:
        raise ValueError("Direction for averaging should be either X, Y or Z!\n")

    tree = ET.parse(potential_file)
    parsed_potential_file = tree.getroot()

    gridticks = [int(n) for n in parsed_potential_file.find("./grid").attrib["gridticks"].split()]
    potential_ = []
    for nr in parsed_potential_file.findall(".//row/row"):
       potential_.append([float(rho) for rho in nr.text.split()])
    potential_ = np.array(list(zip(*[iter(potential_)]*gridticks[1]))).swapaxes(0,2)

    # Remove last points from data (to avoid double counting due to PBC)
    potential = potential_[:-1, :-1, :-1]
    ngridpts = np.array(potential.shape)

    endpointrs = []
    for axis in parsed_potential_file.findall("./grid/axis"):
        endpointrs.append(axis.get("endpointrs"))

    cell = []
    for vec in endpointrs:
        cell.append(vec.split())
    cell = np.array(cell).astype(float)
    latticelength = np.dot(cell, cell.T).diagonal() ** 0.5

    # Perform average
    direction_map = {"X": 0, "Y": 1, "Z": 2}
    idir = direction_map[direction]
    a = (idir + 1) % 3
    b = (idir + 2) % 3

    # At each point, sum over other two indices and scale by number of grid points in the plane
    average = []
    distance = []
    xdiff = latticelength[idir] / ngridpts[idir]
    for ipt in range(ngridpts[idir]):
        if direction == "X":
            average.append(potential[ipt, :, :].sum() / (ngridpts[a] * ngridpts[b]))
        elif direction == "Y":
            average.append(potential[:, ipt, :].sum() / (ngridpts[a] * ngridpts[b]))
        else:
            average.append(potential[:, :, ipt].sum() / (ngridpts[a] * ngridpts[b]))
        distance.append(ipt * xdiff)

    average.append(average[0])
    distance.append(ngridpts[idir] * xdiff)
    planar_average_data = np.vstack((distance, average)).T

    return planar_average_data


def main() -> None:
    parser = ArgumentParser(description="""Extract planar-averaged electrostatic potential in a given direction.""")

    parser.add_argument("--potential_file", "-f",
                        type = Union[str, pathlib.Path],
                        default = ["VCL3D.xml"],
                        nargs = 1,
                        dest = "pot_file",
                        help = "name of file containing electrostatic potential data")

    parser.add_argument("direction",
                        type=str,
                        default=["z"],
                        nargs=1,
                        help="direction along which potential needs to be averaged")

    args = parser.parse_args()

    planar_average_data = execute_planar_average(args.pot_file[0], args.direction[0])

    with open(f"planarAverage_{args.direction[0].lower()}", "w") as f:
        np.savetxt(f, planar_average_data, fmt="%15.8f %15.8f")

if __name__ == "__main__":
    main()
