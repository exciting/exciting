"""Set path to species files in **input.xml** files.

Call as:

```bash
python3 -m excitingscripts.setup.speciespath
```
"""

import os
from argparse import ArgumentParser 
from excitingtools import ExcitingInputXML
from excitingtools.exciting_dict_parsers import input_parser

def set_exciting_root(input_file: str, 
                      output_file: str, 
                      speciespath='.') -> None:
    """Replace the path to the species files in `input_file`
    and write to `output_file`.

    :param input_file: Input file.
    :param output_file: Input file, with '$EXCITINGROOT' replaced with excitingroot.
    :param excitingroot: Environment variable string.
    """

    input_data = input_parser.parse_input_xml(input_file)
    input_data['structure']['species_path'] = os.path.abspath(str(speciespath))
    ExcitingInputXML(**input_data).write(output_file)

def main() -> None:

    parser = ArgumentParser(description="Set the path to the species files in input files.")

    parser.add_argument("--input-file", "-i", 
                        type = str,
                        default = ["input.xml"],
                        nargs = 1,
                        dest = "infile", 
                        help = "name of the input file")

    parser.add_argument("--output-file", "-o",
                    type = str,
                    nargs = 1,
                    default = ["input.xml"],
                    dest = "outfile",
                    help = "name of the output file")

    parser.add_argument("--path", "-p",
                    type = str,
                    nargs = 1,
                    default = ['.'],
                    dest = "path",
                    help = "Path to species files")

    args = parser.parse_args()

    set_exciting_root(args.infile[0], args.outfile[0], speciespath=args.path[0])

if __name__ == "__main__":
    main()