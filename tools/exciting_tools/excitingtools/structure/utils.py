"""Utils for structure submodule."""

from __future__ import annotations

from pathlib import Path
from typing import Dict

from excitingtools import ExcitingStructure
from excitingtools.species import SpeciesFile


def get_species_symbols(
    structure: ExcitingStructure,
    *,
    read_species_files: bool = False,
    species_files: Dict[str, SpeciesFile] | None = None,
) -> list[str]:
    """Get the species symbols for an exciting structure object, if desired, from the species files.

    :param structure: input exciting structure object
    :param read_species_files: if True, read species files to get the species symbols instead of the filenames.
     Note that the speciespath in the structure is used and should be absolute to guarantee that the path is found.
    :param species_files: dictionary of species files, to take the species symbols from.
     Only one of 'read_species_files' and 'species_files' can be set.
    :return: species symbols
    """
    if not read_species_files and not species_files:
        return structure.species
    if read_species_files and species_files:
        raise ValueError("Only one of 'read_species_files' and 'species_files' can be set.")

    if species_files and not set(structure.unique_species).issubset(species_files):
        raise ValueError("Not all species in the structure are present in the passed species files dict.")
    if read_species_files:
        species_path = Path(structure.structure_attributes["speciespath"]).resolve()
        species_files = {sp: SpeciesFile.from_file(species_path / (sp + ".xml")) for sp in structure.unique_species}

    return [species_files[sp].species["chemicalSymbol"] for sp in structure.species]
