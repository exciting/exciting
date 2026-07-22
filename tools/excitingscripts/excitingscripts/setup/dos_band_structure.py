"""Add DOS and band structure element to given input file by getting the band path from the input structure.

Located at `excitingscripts/setup/dos_band_structure.py`.

Call as:

```bash
python3 -m excitingscripts.setup.dos_band_structure
```
"""

import os
import pathlib
from argparse import ArgumentParser
from os.path import join
from typing import Union

from excitingtools import ExcitingInputXML, ExcitingPropertiesInput
from excitingtools.input.bandstructure import get_bandstructure_input_from_exciting_structure


def setup_dos_band_structure(input_file: Union[str, pathlib.Path], root_directory=os.getcwd()) -> None:
    """Add DOS and band structure element to given input file by getting the band path from the input structure.

    :param input_file: Input file.
    :param root_directory: Root directory.
    """
    parsed_input = ExcitingInputXML.from_xml(input_file)
    dos = {"nsmdos": 2, "nwdos": 1000, "winddos": [-0.3, 0.3]}
    band_structure = get_bandstructure_input_from_exciting_structure(
        parsed_input.structure)  # pylint: disable=no-member
    parsed_input.properties = ExcitingPropertiesInput(dos=dos, bandstructure=band_structure)

    parsed_input.write(join(root_directory, "input.xml"))


def main() -> None:
    parser = ArgumentParser(description="""Add DOS and band structure element to given input file by getting the band 
                                        path from the input structure.""")

    parser.add_argument("--input-file", "-i",
                        type=Union[str, pathlib.Path],
                        default=["input.xml"],
                        nargs=1,
                        dest="infile",
                        help="name of the input file")

    parser.add_argument("--root-directory", "-r",
                        type=Union[str, pathlib.Path],
                        default=[os.getcwd()],
                        nargs=1,
                        dest="root_directory",
                        help="root path for folders that are created by this script")

    args = parser.parse_args()

    setup_dos_band_structure(args.infile[0], args.root_directory[0])


if __name__ == "__main__":
    main()
