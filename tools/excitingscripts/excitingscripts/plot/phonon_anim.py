import os
from argparse import ArgumentParser
from pathlib import Path
from typing import Tuple, Union

import numpy as np
from excitingtools import ExcitingInputXML, parse
from scipy.constants import physical_constants, pi
from excitingtools.exciting_dict_parsers.species_parser import parse_species_xml

# Constants
bohr_to_angstrom = 0.52917721092
hartree_to_cm1 = physical_constants["hartree-inverse meter relationship"][0] * 1e-2
two_pi = 2 * pi


def parse_species_data(unique_species: list, excitingroot: Union[str, Path], all_species: list) -> Tuple[list, int]:
    """Parse data for each species to extract mass, atomic number, and count information.

    :param unique_species: List of unique species in the structure.
    :param excitingroot: exciting root directory.
    :param all_species: List of all species present in the structure.
    :return: Species data containing chemical symbol, mass, atomic number and number of atoms per species, as well as
    maximum atom count of any species in the structure.
    """
    species_data = []
    natmax = 0
    for species in unique_species:
        species_file = Path(excitingroot) / "species" / f"{species}.xml"
        parsed_species = parse_species_xml(species_file)
        count = all_species.count(species)
        species_data.append({
            'name': species,
            'mass': parsed_species["species"]["mass"],
            'atomic_number': int(-parsed_species["species"]["z"]),
            'count': count
        })
        natmax = max(natmax, count)
    return species_data, int(natmax)


def generate_visualization_files(supercell_dims: Tuple[int, int, int], scaling: float, nsteps: int,
                                 root_dir: Union[str, Path]) -> None:
    """Generate .axsf and .xyz files for visualizing phonon modes.

    :param supercell_dims: Dimensions of the supercell as (n1, n2, n3).
    :param scaling: Scaling factor for the atomic displacements.
    :param nsteps: Number of steps in the animation sequence.
    :param root_dir: Directory containing the files input.xml and PHONON.OUT.
    """
    input_file = f"{root_dir}/input.xml"
    parsed_input = ExcitingInputXML.from_xml(input_file)

    atom_positions = parsed_input.structure.positions

    scale = getattr(parsed_input.structure.crystal_properties, "scale", 1.0)
    stretch = getattr(parsed_input.structure.crystal_properties, "stretch", [1.0, 1.0, 1.0])
    lattice_vectors = parsed_input.structure.lattice * scale * bohr_to_angstrom * stretch

    cartesian = getattr(parsed_input.structure, "cartesian", False)

    excitingroot = os.getenv("EXCITINGROOT")
    species_data, natmax = parse_species_data(parsed_input.structure.unique_species, excitingroot,
                                              parsed_input.structure.species)

    phonon_data = parse(f"{root_dir}/PHONON.OUT")

    for iqpt, q_point_data in phonon_data.items():
        q_vector = q_point_data['q_vector']
        for mode in q_point_data['modes']:
            output_axsf = Path(root_dir) / f'q{iqpt}_mode{mode["mode_index"]}.axsf'
            output_xyz = Path(root_dir) / f'q{iqpt}_mode{mode["mode_index"]}.xyz'

            with output_axsf.open('w') as axsf_file, output_xyz.open('w') as xyz_file:
                axsf_file.write(f"ANIMSTEPS {nsteps}\n")
                for istep in range(1, nsteps + 1):
                    axsf_file.write(f"ATOMS {istep}\n")
                    xyz_file.write(f"{len(atom_positions) * np.prod(supercell_dims)}\n")
                    xyz_file.write(f"Step {istep} of {nsteps}\n")

                    for n1, n2, n3 in np.ndindex(supercell_dims):
                        pos_shifted = np.zeros((len(species_data), natmax, 3))
                        for atom_data in mode['eigenvector_info']:
                            species_idx = atom_data['species'] - 1
                            atom_idx = atom_data['atom'] - 1
                            polarisation = atom_data['polarisation'] - 1
                            mass = species_data[species_idx]['mass']
                            disp_real = atom_data['eigenvector_component_real'] / np.sqrt(mass) * scaling
                            disp_imag = atom_data['eigenvector_component_imag'] / np.sqrt(mass) * scaling
                            qR = np.dot(q_vector, [n1, n2, n3])

                            R = np.dot([n1, n2, n3], lattice_vectors[polarisation])
                            if cartesian:
                                pos = atom_positions[atom_idx, polarisation] * bohr_to_angstrom
                            else:
                                pos = np.dot(atom_positions[atom_idx], lattice_vectors[polarisation])

                            pos_shifted[species_idx, atom_idx, polarisation] = (
                                    pos + R +
                                    disp_real * np.cos(two_pi * (qR - istep / nsteps)) -
                                    disp_imag * np.sin(two_pi * (qR - istep / nsteps))
                            )

                        for species in species_data:
                            isp = species_data.index(species)
                            for ia in range(species['count']):
                                coords = pos_shifted[isp, ia, :]
                                xyz_file.write(f"{species['name']} {coords[0]} {coords[1]} {coords[2]}\n")
                                axsf_file.write(
                                    f"{species['atomic_number']} {coords[0]} {coords[1]} {coords[2]} 0.0 0.0 0.0\n")
def main():
    parser = ArgumentParser(description="Generate files for phonon mode visualization.")

    parser.add_argument("--root-directory", "-r",
                        default=os.getcwd(),
                        help="Root directory containing input files.")

    parser.add_argument("--supercell", "-s",
                        nargs=3,
                        type=int,
                        default=(6, 6, 6),
                        help="Supercell dimensions n1, n2, n3")

    parser.add_argument("--scaling", "-c",
                        type=float,
                        default=10.0,
                        help="Scaling factor for displacements.")

    parser.add_argument("--nsteps", "-n",
                        type=int,
                        default=20,
                        help="Number of steps for the animation.")

    args = parser.parse_args()
    supercell_dims = tuple(args.supercell)
    generate_visualization_files(supercell_dims, args.scaling, args.nsteps, args.root_directory)

if __name__ == "__main__":
    main()
