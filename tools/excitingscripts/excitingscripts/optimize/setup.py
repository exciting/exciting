"""Python script for generating structures at different volume/strains."""

import shutil
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from excitingscripts.lattice.parameters import run_sgroup
from excitingtools import ExcitingInputXML

sgroup_in = "sgroup.in"
sgroup_out = "sgroup.out"
sgroup_err = "sgroup.err"

space_group_mappings = {
    (1, 2): "Triclinic",
    (3, 15): "Monoclinic",
    (16, 74): "Orthorhombic",
    (75, 88): "Tetragonal II",
    (89, 142): "Tetragonal I",
    (143, 148): "Rhombohedral II",
    (149, 167): "Rhombohedral I",
    (168, 176): "Hexagonal II",
    (177, 194): "Hexagonal I",
    (195, 206): "Cubic II",
    (207, 230): "Cubic I",
}

optimization_map = {
    "VOL": "volume",
    "BOA": "b/a ratio with constant volume",
    "COA": "c/a ratio with constant volume",
    "ALPHA": "alpha angle with constant volume",
    "BETA": "beta angle with constant volume",
    "GAMMA": "gamma angle with constant volume",
}

optimization_options = {
    "Triclinic": ["VOL", "BOA", "COA", "ALPHA", "BETA", "GAMMA"],
    "Monoclinic": ["VOL", "BOA", "COA", "GAMMA"],
    "Orthorhombic": ["VOL", "BOA", "COA"],
    "Tetragonal I": ["VOL", "COA"],
    "Tetragonal II": ["VOL", "COA"],
    "Rhombohedral I": ["VOL", "COA"],
    "Rhombohedral II": ["VOL", "COA"],
    "Hexagonal I": ["VOL", "COA"],
    "Hexagonal II": ["VOL", "COA"],
    "Cubic I": ["VOL"],
    "Cubic II": ["VOL"],
}

# Dictionary where each key represents the optimization parameter, and the value is a string
# representing the deformation matrix.

deformation_matrix_string = {
    "VOL": "\n[ 1+eps  0      0     ]"
    "\n[ 0      1+eps  0     ]"
    "\n[ 0      0      1+eps ]",
    "BOA": "\n[(1+eps)^-.5   0           0          ]"
    "\n[ 0           (1+eps)^+1.  0          ]"
    "\n[ 0            0          (1+eps)^-.5 ]",
    "COA": "\n[(1+eps)^-.5   0           0          ]"
    "\n[ 0           (1+eps)^-.5  0          ]"
    "\n[ 0            0          (1+eps)^+1. ]",
    "ALPHA": "\n[ 1/(1-eps^2)  0           0          ]"
    "\n[ 0            1          eps         ]"
    "\n[ 0           eps          1          ]",
    "BETA": "\n[ 1           0           eps         ]"
    "\n[ 0           1/(1-eps^2)  0          ]"
    "\n[eps          0            1          ]",
    "GAMMA": "\n[ 1          eps           0          ]"
    "\n[eps          1            0          ]"
    "\n[ 0           0            1/(1-eps^2)]",
}

def get_deformation_matrix(eps: float) -> dict:
    """Get deformation matrix for a given strain value depending on the optimization parameter.

     :param eps: Strain value.
     :return: Dictionary containing deformation_matrix for each optimization parameter.
     """
    return {
        "VOL": np.diag([1.0 + eps, 1.0 + eps, 1.0 + eps]),
        "BOA": np.diag([(1.0 + eps) ** (-0.5), 1.0 + eps, (1.0 + eps) ** (-0.5)]),
        "COA": np.diag([(1.0 + eps) ** (-0.5), (1.0 + eps) ** (-0.5), 1.0 + eps]),
        "ALPHA": np.array(
            [[1.0 / (1.0 - eps**2.0), 0.0, 0.0], [0.0, 1.0, eps], [0.0, eps, 1.0]]
        ),
        "BETA": np.array(
            [[1.0, 0.0, eps], [0.0, 1.0 / (1.0 - eps**2.0), 0.0], [eps, 0.0, 1.0]]
        ),
        "GAMMA": np.array(
            [[1.0, eps, 0.0], [eps, 1.0, 0.0], [0.0, 0.0, 1.0 / (1.0 - eps**2.0)]]
        ),
    }

def get_crystal_system(space_group_number: int) -> str:
    """Get crystal system based on the space group number.

    :param space_group_number: Space group number of the crystal structure.
    :return: Crystal system corresponding to the space group number.
    """
    for space_group_range, crystal_system in space_group_mappings.items():
        if space_group_range[0] <= space_group_number <= space_group_range[1]:
            return crystal_system


def check_monoclinic_compatibility(
    base_vectors: np.ndarray, ref_scale: float, stretch: list, threshold_angle: float
) -> None:
    """
    Checks if the given monoclinic structure is compatible with certain geometric criteria.

    :param base_vectors: The base vectors of the crystal structure.
    :param ref_scale: The reference scale for the unit cell.
    :param stretch: The stretch factors along each axis.
    :param threshold_angle: The threshold for determining if the angle is effectively 90 degrees.
    """

    adjusted_vectors = base_vectors * np.array(stretch)[:, np.newaxis] * ref_scale
    scalar_product = np.dot(adjusted_vectors[0], adjusted_vectors[1])

    # Check if the scalar product is below the threshold, indicating near orthogonality
    if np.abs(scalar_product) < threshold_angle:
        raise ValueError(
            """Your MONOCLINIC structure is not compatible\n
        with the OPTIMIZE internal representation, where\n
        the angle GAMMA (between bvec_1 and bvec_2) is the\n
        ONLY non right angle between the crystal basis vectors!\n
        Please, CHECK your input file!"""
        )


def setup_optimize_lattice(
    run_directory: str, max_strain: float, num_dist_str: int, input_file: str, opt_index: int
) -> None:
    """Set up the optimization of the lattice parameters.

    :param max_strain: The maximum physical strain.
    :param num_dist_str: The number of distorted structures.
    :param infile: Name of input file.
    :param opt_index: The index of the optimization parameter.
    """
    if opt_index < 1 or opt_index > 6:
        raise ValueError("Invalid optimization index.")

    if max_strain < 0.0 or max_strain > 1.0:
        raise ValueError("The maximum physical strain is out of range")

    if num_dist_str % 2 == 0:
        num_dist_str += 1

    if num_dist_str < 4 or num_dist_str > 100:
        raise ValueError(
            "The number of distorted structures must be an odd number greater than 4 and less than 100."
        )

    num_structures_per_side = int(
        (num_dist_str - 1) / 2
    )
    threshold_strain_interval = 0.00001
    if max_strain / num_structures_per_side < threshold_strain_interval:
        raise ValueError(
            """The maximum physical strain is too small for the number of distorted structures.\n
             Choose a larger value for maximum physical strain or less number of distorted structures."""
        )

    # Create a temporary directory to store sgroup.in, sgroup.out, and sgroup.err files
    temp_dir = run_directory / "temp"

    run_sgroup(input_file, temp_dir)

    s_out = temp_dir / sgroup_out
    s_err = temp_dir / sgroup_err

    if s_err.exists() and s_err.stat().st_size != 0:
        raise RuntimeError(f"{sgroup_err} file is not empty. Check the file for errors.")

    if not s_out.exists():
        raise FileNotFoundError(f"{sgroup_out} file not found.")

    lines = s_out.read_text().split("\n")
    for line in lines:
        if line.startswith("Number and name of space group:"):
            space_group_number = int(line.split(":")[1].split("(")[0].strip())
            space_group_symbol = line.split(":")[1].split("(")[1].split(")")[0].strip()
            break

    # Remove temporary directory
    shutil.rmtree(temp_dir)

    crystal_system = get_crystal_system(space_group_number)

    if crystal_system is None:
        raise ValueError(f"Incorrect Space-Group Number: {space_group_number}")

    print(
        f"     Number and name of space group: {space_group_number} ({space_group_symbol}) "
    )
    print(f"     {crystal_system} structure in the Laue classification.\n")


    parsed_input = ExcitingInputXML.from_xml(input_file)

    if hasattr(parsed_input.structure.crystal_properties, "scale"):
        ref_scale = parsed_input.structure.crystal_properties.scale
    else:
        ref_scale = 1.0

    if hasattr(parsed_input.structure.crystal_properties, "stretch"):
        stretch = parsed_input.structure.crystal_properties.stretch
    else:
        stretch = [1.0, 1.0, 1.0]

    base_vectors = np.array(parsed_input.structure.lattice)
    volume = abs(np.linalg.det(base_vectors) * ref_scale**3 * np.prod(stretch))

    # Check compatibility for monoclinic structures
    if crystal_system == "Monoclinic":
        check_monoclinic_compatibility(base_vectors, volume, ref_scale, stretch)

    optimization_parameters = optimization_options[crystal_system]

    # Show the available parameters to optimize
    print("     Parameters to optimize:\n")
    for i, parameter in enumerate(optimization_parameters):
        print(f"     {i + 1} ... {optimization_map[parameter]}")
    print()

    if len(optimization_parameters) < opt_index:
        raise ValueError(
            f"""Invalid optimization index for {crystal_system} structure.\n
             Choose from {optimization_parameters}."""
        )

    setup_dir_name = optimization_parameters[opt_index - 1]

    print(
        f"     {opt_index} -> Parameter to optimize is {optimization_map[setup_dir_name]}\n"
    )
    print(f"     Maximum physical strain is {max_strain}\n")
    print(f"     Number of distorted structures is {num_dist_str}\n")

    setup_dir = run_directory / setup_dir_name

    if setup_dir.exists():
        print(
            f"Warning: Directory {setup_dir_name} already exists. Renaming to {setup_dir_name}_old"
        )
        new_dir_name = f"{setup_dir.name}_old"
        new_dir_path = setup_dir.parent / new_dir_name
        shutil.rmtree(new_dir_path, ignore_errors=True)
        setup_dir.rename(new_dir_path)

    setup_dir.mkdir()
    sourcefile = setup_dir / "source.xml"
    shutil.copy(input_file, sourcefile)

    infofile = setup_dir / f"INFO_{setup_dir_name}"
    with open(infofile, "w") as info:
        info.write(f"Space-group number              = {space_group_number}")
        info.write(f"\nStructure type                  = {crystal_system}")
        info.write(f"\nMaximum physical strain         = {max_strain}")
        info.write(f"\nNumber of distorted structures  = {num_dist_str}")

    paramfile = setup_dir / f"{setup_dir_name.lower()}-Parameters"
    with open(paramfile, "w") as param:
        param.write(
            f"{setup_dir_name}, Deformation Matrix = {deformation_matrix_string[setup_dir_name]}\n"
        )

    # Generate the distorted structures
    i = 0
    for s in range(-num_structures_per_side, num_structures_per_side + 1):

        strain = max_strain * s / num_structures_per_side
        if s == 0:
            strain = threshold_strain_interval

        # Calculate the deformation matrix
        transformation_matrix = get_deformation_matrix(strain)[setup_dir_name]

        # Update the base vectors
        new_base_vectors = np.dot(base_vectors, transformation_matrix)
        parsed_input.structure.lattice = new_base_vectors.tolist()

        strain_dir = setup_dir / f"{setup_dir_name.lower()}_{i + 1}"
        strain_dir.mkdir()
        i = i + 1

        parsed_input.write(strain_dir / "input.xml")
        with open(paramfile, "a") as param:
            param.write(
                f"\n{setup_dir_name.lower()}_{i + 1}, Physical strain = {strain}"
            )
            param.write(
                f"\nV1 --=>    {new_base_vectors[0, 0]:22.16f}    {new_base_vectors[0, 1]:22.16f}    {new_base_vectors[0, 2]:22.16f}"
            )
            param.write(
                f"\nV2 --=>    {new_base_vectors[1, 0]:22.16f}    {new_base_vectors[1, 1]:22.16f}    {new_base_vectors[1, 2]:22.16f}"
            )
            param.write(
                f"\nV3 --=>    {new_base_vectors[2, 0]:22.16f}    {new_base_vectors[2, 1]:22.16f}    {new_base_vectors[2, 2]:22.16f}\n"
            )

    with open(paramfile, "a") as param:
        param.write("\n   Distorted parameters: END")


def main() -> None:
    help_string = "Choose from the following numbers to choose the lattice parameters to be optimized:"
    for crystal_system, parameters in optimization_options.items():
        help_string += f'"{crystal_system}":'
        c = 1
        for parameter in parameters:
            help_string += f"{c} -> {optimization_map[parameter]} , "
            c += 1
        help_string += " | "

    parser = ArgumentParser(
        description="Python script for generating structures at different volume/strains."
    )

    parser.add_argument(
        "-d",
        "--directory",
        type=str,
        default=Path.cwd(),
        dest="run_directory",
        help="Directory where exciting runs",
    )

    parser.add_argument(
        "--input-file", "-i",
        type=str,
        default=str(Path.cwd() / "input.xml"),
        help="Path to the exciting input file."
    )

    parser.add_argument(
        "--opt",
        "-optimization_index",
        type=int,
        default=1,
        dest="optimization_index",
        help=help_string,
    )

    parser.add_argument(
        "--maximum-strain",
        "--s-max",
        type=float,
        required=True,
        dest="max_strain",
        help="Maximum physical strain value. Suggested value is between 0.001 and 0.050.",
    )

    parser.add_argument(
        "--number-of-distorted-structures",
        "-n",
        type=int,
        required=True,
        dest="num_dist_str",
        help="The number of the distorted structures [4 < odd number < 100]",
    )

    args = parser.parse_args()
    run_directory= args.run_directory
    input_file = Path(args.input_file)
    opt_index = args.optimization_index
    max_strain = args.max_strain
    num_dist_str = args.num_dist_str

    # Call the setup_optimize_lattice function
    setup_optimize_lattice(
        run_directory=run_directory,
        max_strain=max_strain,
        num_dist_str=num_dist_str,
        input_file=input_file,
        opt_index=opt_index,
    )


if __name__ == "__main__":
    main()
