""" Test the utilities for ase. """
from ase.build import bulk

from excitingtools import ExcitingStructure
from excitingtools.structure.ase_utilities import exciting_structure_to_ase


def test_class_exciting_structure_to_ase():
    ase_atoms = bulk("Si")
    structure = ExcitingStructure(ase_atoms, species_path='./')
    new_ase_atoms = exciting_structure_to_ase(structure)
    assert ase_atoms.wrap() == new_ase_atoms.wrap()
