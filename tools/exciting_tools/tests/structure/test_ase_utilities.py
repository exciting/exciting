"""Test the utilities for ase."""

import pytest

from excitingtools import ExcitingStructure


def test_class_exciting_structure_to_ase():
    ase_build = pytest.importorskip("ase.build")
    ase_utilities = pytest.importorskip("excitingtools.structure.ase_utilities")

    ase_atoms = ase_build.bulk("Si")
    structure = ExcitingStructure(ase_atoms, species_path="./")
    new_ase_atoms = ase_utilities.exciting_structure_to_ase(structure)
    assert ase_atoms.wrap() == new_ase_atoms.wrap()
