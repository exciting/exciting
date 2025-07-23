"""Generate structures with different interlayer distances.

Located at `excitingscripts/setup/interlayer_distance.py`.

Call as:

```bash
python3 -m excitingscripts.setup.interlayer_distance dmin dmax nr_displ dinfty
```
Where <code>dmin</code> and <code>dmax</code> are the minimum and maximum values for the interlayer distance, <code>nr_displ</code> is the number of distances in the interval [<code>dmin</code>, <code>dmin</code>] and <code>dinfty</code> is the interlayer distance at infinity.
"""


import os
import pathlib
from argparse import ArgumentParser
from os.path import join
from typing import Union

from excitingtools import ExcitingInputXML


def setup_interlayer_distance(input_file: Union[str, pathlib.Path], dmin: float, dmax: float, displ_points: int,
                              dinfty=None, root_directory=os.getcwd()) -> None:
    """Create input files for structures with different interlayer distances and save them in corresponding directories.

    :param input_file: Input file.
    :param dmin: Minimum interlayer distance in Bohr.
    :param dmax: Maximum interlayer distance in Bohr.
    :param displ_points: Number of distances in [dmin, dmax].
    :param dinfty: Interlayer distance at infinity in Bohr.
    :param root_directory: Root directory.
    """
    if not (0 < displ_points < 99):
        raise ValueError("Number of displacements is out of range [0-99]!\n")

    setup_infinity = True
    if dinfty is None or dinfty <= dmax:
        setup_infinity = False

    parsed_input = ExcitingInputXML.from_xml(input_file)

    if hasattr(parsed_input.structure.crystal_properties, "scale"):
        scale = parsed_input.structure.crystal_properties.scale
    else:
        scale = 1.0

    if hasattr(parsed_input.structure.crystal_properties, "stretch"):
        stretch = parsed_input.structure.crystal_properties.stretch
    else:
        stretch = [1.0, 1.0, 1.0]

    delta = displ_points - 1
    displ_step = float(dmax - dmin) / delta

    for i in range(displ_points):
        displ = dmin + i * displ_step

        # name and create directory for storing calculations with stretched cells
        rundir = join(root_directory, f"rundir-{i + 1}")
        os.makedirs(rundir, exist_ok=True)

        with open(join(rundir, f"strain-{i + 1}"), "w") as fid:
            fid.write(str(displ))

        parsed_input.structure.lattice[2][2] = 2 * displ / scale / stretch[2]

        parsed_input.write(join(rundir, "input.xml"))

    if setup_infinity:
        displ = dinfty

        rundir = join(root_directory, "rundir-oo")

        os.makedirs(rundir, exist_ok=True)

        with open(join(rundir, "strain-oo"), "w") as fid:
            fid.write(str(displ))

        parsed_input.structure.lattice[2][2] = 2 * displ / scale / stretch[2]

        parsed_input.write(join(rundir, "input.xml"))

def main() -> None:
    parser = ArgumentParser(description="Create input files for structures with different interlayer distances.")

    parser.add_argument("--input-file", "-i",
                        type=Union[str, pathlib.Path],
                        default=["input.xml"],
                        nargs=1,
                        dest="infile",
                        help="name of the input file")

    parser.add_argument("dmin",
                        type=float,
                        nargs=1,
                        help="minimum interlayer distance")

    parser.add_argument("dmax",
                        type=float,
                        nargs=1,
                        help="maximum interlayer distance")

    parser.add_argument("displ_points",
                        type=int,
                        nargs=1,
                        help="number of distances in [dmin, dmax]")

    parser.add_argument("dinfty",
                        type=float,
                        default=[None],
                        nargs=1,
                        help="interlayer distance at infinity")

    parser.add_argument("--root-directory", "-r",
                        default=[os.getcwd()],
                        nargs=1,
                        dest="root_directory",
                        help="root path for folders that are created by this script")

    args = parser.parse_args()

    setup_interlayer_distance(args.infile[0], args.dmin[0], args.dmax[0], args.displ_points[0], args.dinfty[0],
                              args.root_directory[0])


if __name__ == "__main__":
    main()
