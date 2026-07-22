"""Add band structure element to given input file by getting the band path from the input structure."""

import os
import pathlib
from argparse import ArgumentParser
from os.path import join
from typing import Union

from excitingtools import ExcitingInputXML, ExcitingPropertiesInput
from excitingtools.input.bandstructure import get_bandstructure_input_from_exciting_structure


def setup_band_structure(input_file: Union[str, pathlib.Path], 
                         root_directory: str=os.getcwd(), 
                         create_directories: bool=True,
                         overwrite: bool=False) -> None:
    """Add band structure element to given input file by getting the band path from the input structure.

    :param input_file: Input file.
    :param root_directory: Root directory.
    """
    parsed_input = ExcitingInputXML.from_xml(input_file)
    band_structure = get_bandstructure_input_from_exciting_structure(
        parsed_input.structure)  # pylint: disable=no-member
    parsed_input.properties = ExcitingPropertiesInput(bandstructure=band_structure)

    if create_directories:
        os.makedirs(root_directory, exist_ok=True)
    
    out_path = join(root_directory, "input.xml")
    if not os.path.exists(out_path) or overwrite:
        parsed_input.write(out_path)
    else:
        raise FileExistsError(f'File at {out_path} exists.')


def main() -> None:
    parser = ArgumentParser(description="""Add band structure element to given input file by getting the band path from
                                        the input structure.""")

    parser.add_argument("--input-file", "-i",
                        type=str,
                        default=["input.xml"],
                        nargs=1,
                        dest="infile",
                        help="name of the input file")

    parser.add_argument("--root-directory", "-r",
                        type=str,
                        default=[os.getcwd()],
                        nargs=1,
                        dest="root_directory",
                        help="root path for folders that are created by this script")

    parser.add_argument("--make-dirs",
                        dest="mkdirs",
                        action='store_true',
                        help="Make directories, does not raise an exception if they exist")

    parser.add_argument("--no-overwrite",
                        default=False,
                        dest="no_overwrite",
                        action='store_true',
                        help="Don't overwrite existing files.")

    args = parser.parse_args()

    try:
        setup_band_structure(args.infile[0], 
                             args.root_directory[0], 
                             create_directories=args.mkdirs, 
                             overwrite=not args.no_overwrite)
    except FileExistsError as e:
        print(f"Can't write file: {e} Use --overwrite flag to overwrite.")

if __name__ == "__main__":
    main()
