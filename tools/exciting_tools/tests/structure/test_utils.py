"""Test utils.py"""

import copy
from pathlib import Path

import pytest

from excitingtools import ExcitingStructure
from excitingtools.structure.utils import get_species_symbols


def test_get_species_symbols_no_species_files(structure_H2He: ExcitingStructure) -> None:
    assert get_species_symbols(structure_H2He, read_species_files=False) == ["H_core", "H", "He"]


def test_get_species_symbols_both_specified(structure_H2He: ExcitingStructure, species_files_h_hcore_he) -> None:
    with pytest.raises(ValueError, match="Only one of 'read_species_files' and 'species_files' can be set."):
        get_species_symbols(structure_H2He, read_species_files=True, species_files=species_files_h_hcore_he)


def test_get_species_symbols_from_species_files_key_missing(
    structure_H2He: ExcitingStructure, species_files_h_hcore_he
) -> None:
    new_species_files = copy.deepcopy(species_files_h_hcore_he)
    new_species_files.pop("H")
    with pytest.raises(
        ValueError, match="Not all species in the structure are present in the passed species files dict."
    ):
        get_species_symbols(structure_H2He, read_species_files=False, species_files=new_species_files)


def test_get_species_symbols_from_species_files(structure_H2He: ExcitingStructure, species_files_h_hcore_he) -> None:
    assert get_species_symbols(structure_H2He, read_species_files=False, species_files=species_files_h_hcore_he) == [
        "H",
        "H",
        "He",
    ]


def test_get_species_symbols_from_species_path(
    structure_H2He: ExcitingStructure, species_files_h_hcore_he, tmp_path: Path
) -> None:
    new_structure = copy.deepcopy(structure_H2He)
    mock_species_path = tmp_path / "species_path"
    mock_species_path.mkdir()
    new_structure.structure_attributes["speciespath"] = str(mock_species_path)

    for key, value in species_files_h_hcore_he.items():
        value.write(mock_species_path / f"{key}.xml")

    assert get_species_symbols(new_structure, read_species_files=True) == ["H", "H", "He"]
