"""Utilities for interaction with the ASE library."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Union

import numpy as np
from ase import Atoms

from excitingtools import ExcitingStructure
from excitingtools.constants.units import bohr_to_angstrom
from excitingtools.exciting_obj_parsers.optimization_trajectory import parse_optimization_trajectory
from excitingtools.species import SpeciesFile
from excitingtools.structure.utils import get_species_symbols


def exciting_structure_to_ase(
    structure: ExcitingStructure,
    *,
    read_species_files: bool = False,
    species_files: Dict[str, SpeciesFile] | None = None,
) -> Atoms:
    """Function to extract the physical structure from an exciting structure object
    and transforms it into an ase.atoms.Atoms object.

    :param structure: input exciting structure object
    :param read_species_files: if True, read species files to get the species symbols instead of the filenames.
     Note that the speciespath in the structure is used and should be absolute to guarantee that the path is found.
    :param species_files: dictionary of species files, to take the species symbols from.
     Only one of 'read_species_files' and 'species_files' can be set.
    :returns: ASE Atoms object
    """
    species = get_species_symbols(structure, read_species_files=read_species_files, species_files=species_files)

    lattice = structure.get_lattice(convert_to_angstrom=True)
    if structure.structure_attributes.get("cartesian"):
        positions = np.asarray(structure.positions) * bohr_to_angstrom
        return Atoms(symbols=species, positions=positions, cell=lattice, pbc=True)
    return Atoms(symbols=species, scaled_positions=structure.positions, cell=lattice, pbc=True)


def parse_optimization_trajectory_ase(info_out_file: Union[str, Path]) -> List[Atoms]:
    """Parse the optimization trajectory into a list of ASE Atoms objects.

    :param info_out_file: path to the INFO.OUT file
    :return: list of structures as ASE Atoms objects
    """
    exciting_structs = parse_optimization_trajectory(info_out_file)
    return [exciting_structure_to_ase(exciting_struct) for exciting_struct in exciting_structs]
