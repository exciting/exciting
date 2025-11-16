"""Utilities to interact with the pymatgen library."""

from __future__ import annotations

from typing import Dict

import numpy as np
from pymatgen.core import Structure

from excitingtools import ExcitingStructure
from excitingtools.constants.units import angstrom_to_bohr, bohr_to_angstrom
from excitingtools.species import SpeciesFile
from excitingtools.structure.utils import get_species_symbols


def exciting_structure_to_pymatgen(
    structure: ExcitingStructure,
    *,
    read_species_files: bool = False,
    species_files: Dict[str, SpeciesFile] | None = None,
) -> Structure:
    """Function to extract the physical structure from an exciting structure object
    and transforms it into a pymatgen.core.structure.Structure object.

    :param structure: input exciting structure object
    :param read_species_files: if True, read species files to get the species symbols instead of the filenames.
     Note that the speciespath in the structure is used and should be absolute to guarantee that the path is found.
    :param species_files: dictionary of species files, to take the species symbols from.
     Only one of 'read_species_files' and 'species_files' can be set.
    :returns: pymatgen Structure object
    """
    species = get_species_symbols(structure, read_species_files=read_species_files, species_files=species_files)

    lattice = structure.get_lattice(convert_to_angstrom=True)
    cartesian = structure.structure_attributes.get("cartesian", False)
    positions = np.asarray(structure.positions) * bohr_to_angstrom if cartesian else structure.positions

    return Structure(lattice=lattice, species=species, coords=positions, coords_are_cartesian=cartesian)


def pymatgen_to_exciting_structure(structure: Structure) -> ExcitingStructure:
    """Initialise lattice, species and positions from a pymatgen Structure Object.
    Note: pymatgen works in Angstrom, whereas exciting expects atomic units

    :param structure: pymatgen Structure object.
    :return exciting structure object
    """
    lattice = structure.lattice.matrix * angstrom_to_bohr
    species = [x.symbol.capitalize() for x in structure.species]
    positions = structure.frac_coords
    atoms = [{"species": atom, "position": positions[i]} for i, atom in enumerate(species)]
    return ExcitingStructure(atoms, lattice)
