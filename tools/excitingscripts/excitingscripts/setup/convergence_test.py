"""Generate input files with different values of the main computational parameters.

Call as:

```bash
python3 -m excitingscripts.setup.convergence_test k_i k_f rgkmax_i rgkmax_f
```
Where <code>k_i</code> and <code>k_f</code> are the initial and final k-values for defining the <code><span style="color:green">groundstate</span></code> attribute  <code><span style="color:MediumBlue">ngridk</span></code>, and <code>rgkmax_i</code> and <code>rgkmax_f</code> the  initial and final values for the  <code><span style="color:green">groundstate</span></code> attribute  <code><span style="color:MediumBlue">rgkmax</span></code>."""

import os
import pathlib
from argparse import ArgumentParser
from os.path import join
from typing import Union
import numpy as np

from excitingtools import ExcitingInputXML


def setup_convergence_test(input_file: Union[str, pathlib.Path], 
                           k_initial: int, 
                           k_final: int, 
                           delta_k: int,
                           rgkmax_initial: int,
                           rgkmax_final: int, 
                           delta_rgkmax: int,
                           root_directory=os.getcwd()) -> None:
    """Create input files with varying values for the groundstate attributes ngridk and rgkmax and save them in
    corresponding directories.

        :param input_file: Input file.
        :param k_initial: Initial k-value for defining the groundstate attribute ngridk.
        :param k_final: Final k-value for defining the groundstate attribute ngridk.
        :param rgkmax_initial: Initial value for the groundstate attribute rgkmax.
        :param rgkmax_final: Final value for the groundstate attribute rgkmax.
        :param root_directory: Root directory.
    """

    parsed_input = ExcitingInputXML.from_xml(input_file)

    for k in range(k_initial, k_final + 1, delta_k):
        parsed_input.groundstate.ngridk = [k, k, k]

        for rgkmax in np.arange(rgkmax_initial, rgkmax_final + (1 * delta_rgkmax), delta_rgkmax):
            parsed_input.groundstate.rgkmax = rgkmax

            out_path = join(root_directory, f"{k}_{rgkmax}")

            os.makedirs(out_path, exist_ok=True)

            parsed_input.write(join(out_path, "input.xml"))

def main() -> None:
    parser = ArgumentParser(description="""Create input files with varying values for the groundstate attributes ngridk
                            and rgkmax and save them in corresponding directories.""")

    parser.add_argument("--input-file", "-i",
                        type=str,
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
    
    parser.add_argument('--delta-k', "-dk",
                        type=int,
                        default=2,
                        
                        help="change of k-grid value")

    parser.add_argument("rgkmax_initial",
                        type=int,
                        nargs=1,
                        help="initial k-value for rgkmax")

    parser.add_argument("rgkmax_final",
                        type=int,
                        nargs=1,
                        help="final k-value for rgkmax")

    parser.add_argument('--delta-rgkmax', "-dr",
                        type=float,
                        default=1,
                        help="change of k-grid value")

    parser.add_argument("--root-directory", "-r",
                        type=str,
                        default=[os.getcwd()],
                        nargs=1,
                        dest="root_directory",
                        help="root path for folders that are created by this script")

    args = parser.parse_args()

    setup_convergence_test(args.infile[0], 
                           args.k_initial[0], 
                           args.k_final[0],
                           args.delta_k, 
                           args.rgkmax_initial[0],
                           args.rgkmax_final[0],
                           args.delta_rgkmax, 
                           args.root_directory[0])


if __name__ == "__main__":
    main()
