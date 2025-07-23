"""Convert **xml** files to **xsf**.

Located at `excitingscripts/convert_xml2xsf.py`.

Call as:

```bash
python3 -m excitingscripts.convert_xml2xsf -f file -d dimension
```
Where <code>file</code> is the **xml** file to be converted to **xsf** and <code>dimension</code> is the dimension of <code><span style="color:green">plot</span></code> sub-element in the <code><span style="color:green">properties</span></code> element for a given exciting calculation.
"""

import os
import subprocess
from argparse import ArgumentParser
from os.path import join


def convert_xml2xsf(file_to_convert: str, dimension: str, excitingroot=os.getenv("EXCITINGROOT")) -> None:
    """Convert a given XML file to XSF.

    :param file_to_convert: XML File to convert.
    :param dimension: Dimension of "plot" sub-element in "properties" element for the given exciting calculation.
    :param excitingroot: Environment variable string.
    """
    excitingvisual = join(excitingroot, "xml/visualizationtemplates")
    excitingconvert = join(excitingroot, "xml/inputfileconverter")

    dimension = dimension.lower()

    if file_to_convert == "input.xml":
        xsltproc_command = f"xsltproc {excitingconvert}/xmlinput2xsf.xsl input.xml > input.xsf"

    else:
        xsltproc_command = f"""
        xsltproc {excitingvisual}/plot{dimension}2xsf.xsl {file_to_convert} > {file_to_convert.replace(".xml", ".xsf")}"""

    process = subprocess.Popen(xsltproc_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    output, errors = process.communicate()

    if errors != "":
        raise RuntimeError(
            f"xsltproc failed with following error: {errors}")

def main() -> None:

    parser = ArgumentParser(description="Convert a given XML file to XSF.")

    parser.add_argument("--convert-file", "-f",
                        type = str,
                        default = ["input.xml"],
                        nargs = 1,
                        dest = "convertfile",
                        help = "name of file to convert")

    parser.add_argument("--dimension", "-d",
                        type=str,
                        default=["3d"],
                        nargs=1,
                        dest="dimension",
                        help="dimension")

    args = parser.parse_args()

    convert_xml2xsf(args.convertfile[0], args.dimension[0])

if __name__ == "__main__":
    main()