"""Replace placeholder "$EXCITINGROOT" in **input.xml** files by actual path.

Located at `excitingscripts/setup/excitingroot.py`.

Call as:

```bash
python3 -m excitingscripts.setup.excitingroot
```
"""

import os
from argparse import ArgumentParser 

def set_exciting_root(input_file: str, 
                      output_file: str, 
                      excitingroot=os.getenv("EXCITINGROOT")) -> None:
    """Replace all instances of the string '$EXCITINGROOT' in the file given by `input_file`
    and write to `output_file`.

    :param input_file: Input file.
    :param output_file: Input file, with '$EXCITINGROOT' replaced with excitingroot.
    :param excitingroot: Environment variable string.
    """
    if not excitingroot:
        raise ValueError(
            "EXCITINGROOT is not defined as an environment variable in the shell.\n"
            "If using bash please type: `export EXCITINGROOT=<path-to-exciting_smp>`")
    
    with open(input_file, "r") as f:
        exciting_input = f.read()
    
    exciting_input = exciting_input.replace("$EXCITINGROOT", excitingroot)

    with open(output_file, "w") as f:
        f.write(exciting_input)


def main() -> None:

    parser = ArgumentParser(description="Set the evironmental valriable `$EXCITINGROOT` in input files.")

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
                    default = [os.getenv("EXCITINGROOT")],
                    dest = "path",
                    help = "Path (excitingroot)")

    args = parser.parse_args()

    set_exciting_root(args.infile[0], args.outfile[0], excitingroot=args.path[0])

if __name__ == "__main__":
    main()