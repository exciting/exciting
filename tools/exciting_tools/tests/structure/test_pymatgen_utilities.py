"""Tests for pymatgen utilities."""

import numpy as np
import pytest

from excitingtools import ExcitingStructure


@pytest.fixture
def pymatgen_atoms_H2O():
    """
    H20 molecule in a big box (angstrom), in pymatgen Structure()
    Converts a List[dict] to pymatgen.core.structure.Structure.
    """
    pymatgen_struct = pytest.importorskip("pymatgen.core.structure")
    cubic_cell = np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]])
    atoms = [
        {"species": "H", "position": [0.00000, 0.75545, -0.47116]},
        {"species": "O", "position": [0.00000, 0.00000, 0.11779]},
        {"species": "H", "position": [0.00000, 0.75545, -0.47116]},
    ]

    symbols = [atom["species"] for atom in atoms]
    positions = [atom["position"] for atom in atoms]
    return pymatgen_struct.Structure(lattice=cubic_cell, species=symbols, coords=positions, coords_are_cartesian=True)


def test_class_exciting_structure_pymatgen(pymatgen_atoms_H2O):
    """
    Test the pymatgen Structure object gets used correctly by the ExcitingStructure constructor.
    """
    pymatgen_conversion = pytest.importorskip("excitingtools.structure.pymatgen_utilities")
    structure = pymatgen_conversion.pymatgen_to_exciting_structure(pymatgen_atoms_H2O)

    assert structure.species == ["H", "O", "H"]
    assert np.allclose(
        structure.lattice,
        [[18.897261246257703, 0.0, 0.0], [0.0, 18.897261246257703, 0.0], [0.0, 0.0, 18.897261246257703]],
    ), "Expect lattice vectors to match input values"

    assert np.allclose(structure.positions, pymatgen_atoms_H2O.frac_coords), "Expect positions to match input values."

    # This just confirms the XML tree is built, not that it is correct.
    xml_structure = structure.to_xml()
    assert list(xml_structure.keys()) == ["speciespath"], "Only expect speciespath in structure xml keys"


def test_class_exciting_structure_to_pymatgen(pymatgen_atoms_H2O):
    pymatgen_conversion = pytest.importorskip("excitingtools.structure.pymatgen_utilities")
    structure = pymatgen_conversion.pymatgen_to_exciting_structure(pymatgen_atoms_H2O)
    new_pymatgen_atoms = pymatgen_conversion.exciting_structure_to_pymatgen(structure)
    assert pymatgen_atoms_H2O == new_pymatgen_atoms


def test_convert_exciting_to_pymatgen_failure(structure_H2He: ExcitingStructure) -> None:
    pymatgen_conversion = pytest.importorskip("excitingtools.structure.pymatgen_utilities")
    with pytest.raises(IndexError):
        # a bit weird that pymatgen doesn't give a nice error here
        # maybe they expect some more properties, but the origin here should be that "H_core" is not in the PSE
        pymatgen_conversion.exciting_structure_to_pymatgen(structure_H2He)


def test_convert_exciting_to_pymatgen_read_species_files(
    structure_H2He: ExcitingStructure, species_files_h_hcore_he
) -> None:
    pymatgen_conversion = pytest.importorskip("excitingtools.structure.pymatgen_utilities")
    pymatgen_struct = pymatgen_conversion.exciting_structure_to_pymatgen(
        structure_H2He, species_files=species_files_h_hcore_he
    )
    assert [x.symbol for x in pymatgen_struct.species] == ["H", "H", "He"]


def test_get_n_irreducible_mesh_points():
    pymatgen_conversion = pytest.importorskip("excitingtools.structure.pymatgen_utilities")
    cubic_lattice = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
    cu_atom = [{"species": "Cu", "position": [0, 0, 0]}]
    struct = ExcitingStructure(cu_atom, cubic_lattice, "./")
    n_irreducible_mesh_points = pymatgen_conversion.get_num_irreducible_k_points(struct, (10, 10, 10))
    assert n_irreducible_mesh_points == 56
