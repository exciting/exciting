from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import numpy as np
import os


def xml2xsf(input_file_name: str = "input.xml", output_file_name: str = "structure.xsf", root_path: str = "."):
    # Parse the input.xml file
    tree = ET.parse(os.path.join(root_path, input_file_name))
    root = tree.getroot()

    # Constant parameter Bohr to Angstrom
    bohr2angstrom = 0.529177249

    # Extract crystal lattice vectors
    crystal = root.find('.//crystal')
    scale = float(crystal.get('scale', 1.0))  # Default to 1.0 if scale is not specified
    lattice = []
    for vec in crystal.findall('basevect'):
        vec = [float(x) * scale * bohr2angstrom for x in vec.text.split()]
        lattice.append(vec)
    lattice = np.array(lattice)

    # Extract atomic positions and symbols
    positions = []
    symbols = []
    cartesian = root.find('.//structure').get('cartesian', 'false').lower() == 'true'

    for species in root.findall('.//species'):
        species_file = species.get('speciesfile')
        if species_file:  # Skip empty species tags
            symbol = species_file.split('.')[0]  # Extract element (e.g., 'Li' from 'Li.xml')
            # print(symbol)
            for atom in species.findall('atom'):
                coord = [float(x) for x in atom.get('coord').split()]
                if cartesian:
                    coord = np.array(coord) * scale * bohr2angstrom  # Convert to Cartesian 
                else:
                    coord = np.array(coord)  # Fractional coordinates
                positions.append(coord)
                symbols.append(symbol)

    # Write CRYSTAL (XSF) format manually 
    with open(os.path.join(root_path, output_file_name), "w") as f:

        # Write CRYSTAL section
        f.write("CRYSTAL\n")
        f.write("PRIMVEC\n")
        for vec in lattice:
            f.write(f"{vec[0]:16.10f} {vec[1]:16.10f} {vec[2]:16.10f}\n")

        f.write("PRIMCOORD\n")
        f.write(f"{len(symbols)} 1\n")      # 1 = no periodic image

        # Write atoms
        for s, pos in zip(symbols, positions):
            f.write(f"{s:2s}  {pos[0]:16.10f} {pos[1]:16.10f} {pos[2]:16.10f}\n")

def main() -> None:
    parser = ArgumentParser(description="""Convert input.xml to structure.xsf.""")

    parser.add_argument("--root-directory", "-r",
                        default=os.getcwd(),
                        type=str,
                        dest="root_directory",
                        help="root path for the input and output")

    parser.add_argument("--input", "-i",
                        default="input.xml",
                        type=str,
                        dest="infile",
                        help="root path for the input and output")

    parser.add_argument("--output", "-o",
                        default="structure.xsf",
                        type=str,
                        dest="outfile",
                        help="root path for the input and output")


    args = parser.parse_args()

    xml2xsf(input_file_name=args.infile, 
            output_file_name=args.outfile, 
            root_path=args.root_directory)


if __name__ == "__main__":
    main()
