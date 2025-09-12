"""Reusable structure and species files definitions."""

from typing import Dict

import pytest

from excitingtools import ExcitingStructure
from excitingtools.species import SpeciesFile


@pytest.fixture
def structure_H2He() -> ExcitingStructure:
    """Structure object initialised with a mock crystal, using mandatory arguments only."""
    cubic_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    arbitrary_atoms = [
        {"species": "H_core", "position": [0, 0, 0]},
        {"species": "H", "position": [1, 0, 0]},
        {"species": "He", "position": [2, 0, 0]},
    ]
    return ExcitingStructure(arbitrary_atoms, cubic_lattice, "./")


@pytest.fixture
def species_files_h_hcore_he() -> Dict[str, SpeciesFile]:
    hydrogen = SpeciesFile(
        species={"chemicalSymbol": "H", "mass": 1.1, "name": "hydrogen", "z": -1.0},
        muffin_tin={},
        atomic_states=[],
        basis={"default": []},
    )
    hydrogen_core = SpeciesFile(
        species={"chemicalSymbol": "H", "mass": 1.1, "name": "hydrogen", "z": -1.0},
        muffin_tin={},
        atomic_states=[],
        basis={"default": []},
    )
    helium = SpeciesFile(
        species={"chemicalSymbol": "He", "mass": 4.2, "name": "helium", "z": -2.0},
        muffin_tin={},
        atomic_states=[],
        basis={"default": []},
    )
    return {"H": hydrogen, "H_core": hydrogen_core, "He": helium}
