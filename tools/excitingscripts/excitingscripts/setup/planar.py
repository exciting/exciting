"""Python script for setting calculations for two-dimensional materials."""

import os
import numpy as np
from argparse import ArgumentParser
from os.path import join
from excitingtools import ExcitingInputXML


def setup_planar(workdir=os.getcwd()) -> None:
    """Creates the planar files

    :param workdir: Working directory.
    """
    source_file = join(workdir, "source.xml")

    if not os.path.exists(source_file):
        raise FileNotFoundError(source_file)

    parsed_input = ExcitingInputXML.from_xml(source_file)

    # Extract crystal properties
    if hasattr(parsed_input.structure.crystal_properties, "scale"):
        ref_scale = parsed_input.structure.crystal_properties.scale
    else:
        ref_scale = 1.0

    base_vectors = np.array(parsed_input.structure.lattice)

    with open(join(workdir, "planar"), 'w') as f:
        f.write(f"{ref_scale} {base_vectors[-1, -1]}")


def main() -> None:
    parser = ArgumentParser(description="Python script for setting calculations for two-dimensional materials.")

    parser.add_argument("--root-directory", "-r",
                        type=str,
                        default=[os.getcwd()],
                        nargs=1,
                        dest="workdir",
                        help="path for folder which contains source.xml")

    args = parser.parse_args()

    workdir = args.workdir[0]

    # Set up the planar file
    setup_planar(workdir)


if __name__ == "__main__":
    main()
