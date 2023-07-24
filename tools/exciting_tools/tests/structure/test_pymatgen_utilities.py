""" Tests for pymatgen utilities. """
import numpy as np
import pytest


@pytest.fixture
def pymatgen_atoms_H2O():
    """
    H20 molecule in a big box (angstrom), in pymatgen Structure()
    Converts a List[dict] to pymatgen.core.structure.Structure.
    """
    pymatgen_struct = pytest.importorskip("pymatgen.core.structure")
    cubic_cell = np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]])
    atoms = [{'species': 'H', 'position': [0.00000, 0.75545, -0.47116]},
             {'species': 'O', 'position': [0.00000, 0.00000, 0.11779]},
             {'species': 'H', 'position': [0.00000, 0.75545, -0.47116]}]

    symbols = [atom['species'] for atom in atoms]
    positions = [atom['position'] for atom in atoms]
    return pymatgen_struct.Structure(lattice=cubic_cell, species=symbols, coords=positions, coords_are_cartesian=True)


def test_class_exciting_structure_pymatgen(pymatgen_atoms_H2O):
    """
    Test the pymatgen Structure object gets used correctly by the ExcitingStructure constructor.
    """
    pymatgen_conversion = pytest.importorskip("excitingtools.structure.pymatgen_utilities")
    structure = pymatgen_conversion.pymatgen_to_exciting_structure(pymatgen_atoms_H2O)

    assert structure.species == ["H", "O", "H"]
    assert np.allclose(structure.lattice,
                       [[18.897261246257703, 0.0, 0.0],
                        [0.0, 18.897261246257703, 0.0],
                        [0.0, 0.0, 18.897261246257703]]), \
        'Expect lattice vectors to match input values'

    assert np.allclose(structure.positions, pymatgen_atoms_H2O.frac_coords), 'Expect positions to match input values.'

    # This just confirms the XML tree is built, not that it is correct.
    xml_structure = structure.to_xml()
    assert list(xml_structure.keys()) == ['speciespath'], 'Only expect speciespath in structure xml keys'


def test_class_exciting_structure_to_pymatgen(pymatgen_atoms_H2O):
    pymatgen_conversion = pytest.importorskip("excitingtools.structure.pymatgen_utilities")
    structure = pymatgen_conversion.pymatgen_to_exciting_structure(pymatgen_atoms_H2O)
    new_pymatgen_atoms = pymatgen_conversion.exciting_structure_to_pymatgen(structure)
    assert pymatgen_atoms_H2O == new_pymatgen_atoms
