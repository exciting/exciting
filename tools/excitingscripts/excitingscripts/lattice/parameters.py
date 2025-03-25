import subprocess
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from excitingtools import ExcitingInputXML

sgroup_in = "sgroup.in"
sgroup_out = "sgroup.out"
sgroup_err = "sgroup.err"


def convert_exciting_to_sgroup(input_file: str, output_file: str) -> None:
    """
    Convert an exciting input file  into an sgroup input file.

    :param input_file: Path to the exciting input.xml file.
    :param output_file: Path to the sgroup input file.
    """
    input_path = Path(input_file)
    output_path = Path(output_file)

    if not input_path.exists():
        raise FileNotFoundError(f"ERROR: exciting input file '{input_file}' not found!")

    parsed_input = ExcitingInputXML.from_xml(input_file)


    scale = getattr(parsed_input.structure.crystal_properties, "scale", 1.0)
    lattice_vectors = parsed_input.structure.lattice * scale

    a, b, c = np.linalg.norm(lattice_vectors, axis=1)
    alpha = np.degrees(np.arccos(np.dot(lattice_vectors[1], lattice_vectors[2]) / (b * c)))
    beta = np.degrees(np.arccos(np.dot(lattice_vectors[0], lattice_vectors[2]) / (a * c)))
    gamma = np.degrees(np.arccos(np.dot(lattice_vectors[0], lattice_vectors[1]) / (a * b)))

    with open(output_path, "w") as output_file:
        output_file.write("P\n")
        output_file.write(f"{a:12.10f} {b:12.10f} {c:12.10f} {alpha:12.10f} {beta:12.10f} {gamma:12.10f}\n\n")

        species = parsed_input.structure.species
        atom_positions = parsed_input.structure.positions

        atom_count = len(species)
        output_file.write(f"{atom_count}\n")

        for atom_index in range(atom_count):
            atom_coords = np.array(atom_positions[atom_index])

            cartesian = getattr(parsed_input.structure, "cartesian", False)
            if cartesian:
                atom_coords = np.linalg.inv(np.transpose(lattice_vectors)) @ atom_coords

            output_file.write(f"{atom_coords[0]:15.9f} {atom_coords[1]:15.9f} {atom_coords[2]:15.9f}\n")
            output_file.write(f"{species[atom_index]}.xml\n")


def run_sgroup(input_file: str, output_dir: str) -> None:
    """Run the sgroup program to extract the space group information from the input file.

    :param input_file: Path to the exciting input.xml file.
    :param output_dir: Path to the output directory.
    """
    input_path = Path(input_file)
    output_dir = Path(output_dir)

    if not input_path.exists():
        raise FileNotFoundError(f"File {input_file} not found.")

    output_dir.mkdir(parents=True, exist_ok=True)

    s_in = output_dir / sgroup_in
    s_out = output_dir / sgroup_out
    s_err = output_dir / sgroup_err

    convert_exciting_to_sgroup(str(input_path), str(s_in))

    with open(s_err, "w") as err_file:
        subprocess.run(["sgroup", str(s_in), str(s_out)], stderr=err_file, check=True)


def get_parameters(directory: str):
    """
    Extract lattice symmetry and parameters from the sgroup.out file.

    :param directory: Directory containing the sgroup.out file.
    :return: Dictionary with lattice symmetry and parameters.
    """
    directory = Path(directory)
    s_out = directory / sgroup_out
    s_err = directory / sgroup_err

    if s_err.exists() and s_err.stat().st_size != 0:
        raise RuntimeError(f"{sgroup_err} file is not empty. Check the file for errors.")

    if not s_out.exists():
        raise FileNotFoundError(f"{sgroup_out} file not found.")

    lines = s_out.read_text().splitlines()
    sym = lines[0].strip()
    abc = list(map(float, lines[3].strip().split()))
    ABG = list(map(float, lines[5].strip().split()))

    return {
        "sym": sym,
        "a": abc[0],
        "b": abc[1],
        "c": abc[2],
        "alpha": ABG[0],
        "beta": ABG[1],
        "gamma": ABG[2],
    }


def print_parameters(parameters: dict) -> None:
    """Print the space group parameters in a formatted way.

    :param parameters: Dictionary with lattice symmetry and parameters.
    """
    print(f"     {parameters['sym']}")
    print("     a          b          c          alpha      beta       gamma")
    print(
        f"{parameters['a']:12.5f} {parameters['b']:10.5f} {parameters['c']:10.5f} "
        f"{parameters['alpha']:10.5f} {parameters['beta']:10.5f} {parameters['gamma']:10.5f}"
    )


def main() -> None:
    parser = ArgumentParser(
        description="Convert an exciting input file to an sgroup input file and extract lattice symmetry and parameters."
    )

    parser.add_argument(
        "--input-file", "-i",
        type=str,
        default=str(Path.cwd() / "input.xml"),
        help="Path to the exciting input file."
    )

    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        default=str(Path.cwd()),
        help="Directory to save output files."
    )

    args = parser.parse_args()

    run_sgroup(args.input_file, args.output_dir)
    parameters = get_parameters(args.output_dir)
    print_parameters(parameters)


if __name__ == "__main__":
    main()
