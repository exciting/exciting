import os
import pathlib
from argparse import ArgumentParser
from os.path import join
from typing import Union

from excitingtools import ExcitingInputXML


def setup_convergence_test(input_file: Union[str, pathlib.Path], k_initial: int, k_final: int, rgkmax_initial: int,
                           rgkmax_final: int, root_directory=os.getcwd()) -> None:
    """Create input files with varying values for the groundstate attributes ngridk and rgkmax and save them in
    corresponding directories.

        :param input_file: Input file.
        :param k_initial: Initial k-value for defining the groundstate attribute ngridk.
        :param k_final: Final k-value for defining the groundstate attribute ngridk.
        :param rgkmax_initial: Initial value for the groundstate attribute rgkmax.
        :param rgkmax_final: Final value for the groundstate attribute rgkmax.
        :param root_directory: Root directory.
    """
    dk = 2

    parsed_input = ExcitingInputXML.from_xml(input_file)

    for k in range(k_initial, k_final + 1, dk):
        parsed_input.groundstate.ngridk = [k, k, k]

        for rgkmax in range(rgkmax_initial, rgkmax_final + 1):
            parsed_input.groundstate.rgkmax = rgkmax

            out_path = join(root_directory, f"{k}_{rgkmax}")

            os.makedirs(out_path, exist_ok=True)

            parsed_input.write(join(out_path, "input.xml"))

def main() -> None:
    parser = ArgumentParser(description="""Create input files with varying values for the groundstate attributes ngridk
                            and rgkmax and save them in corresponding directories.""")

    parser.add_argument("--input-file", "-i",
                        type=Union[str, pathlib.Path],
                        default=["input.xml"],
                        nargs=1,
                        dest="infile",
                        help="name of the input file")

    parser.add_argument("k_initial",
                        type=int,
                        nargs=1,
                        help="initial k-value for defining ngridk")

    parser.add_argument("k_final",
                        type=int,
                        nargs=1,
                        help="final k-value for defining ngridk")

    parser.add_argument("rgkmax_initial",
                        type=int,
                        nargs=1,
                        help="initial k-value for rgkmax")

    parser.add_argument("rgkmax_final",
                        type=int,
                        nargs=1,
                        help="final k-value for rgkmax")

    parser.add_argument("--root-directory", "-r",
                        type=Union[str, pathlib.Path],
                        default=[os.getcwd()],
                        nargs=1,
                        dest="root_directory",
                        help="root path for folders that are created by this script")

    args = parser.parse_args()

    setup_convergence_test(args.infile[0], args.k_initial[0], args.k_final[0], args.rgkmax_initial[0],
                           args.rgkmax_final[0], args.root_directory[0])


if __name__ == "__main__":
    main()
