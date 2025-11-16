"""Utilities for interaction with the ASE library."""

from __future__ import annotations

from typing import Dict

import numpy as np
from ase import Atoms

from excitingtools import ExcitingStructure
from excitingtools.constants.units import bohr_to_angstrom
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
