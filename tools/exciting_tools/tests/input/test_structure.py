"""Test ExcitingStructure, python API that generates exciting's structure XML.

NOTE:
All attribute tests should assert on the XML tree content,s as the attribute
order is not preserved by the ElementTree.tostring method. Elements appear to
be fine.

For example:

 `gs_xml_string = xml.etree.ElementTree.tostring(
     gs_xml, encoding='unicode', method='xml')`

may return:

'<groundstate ngridk="8 8 8" rgkmax="8.6"> </groundstate>'
or
'<groundstate rgkmax="8.6" ngridk="8 8 8"> </groundstate>'

"""
import sys
from typing import List
import py
import pytest
import numpy as np

try:
    import ase
except ImportError:
    pass

from excitingtools.input.structure import ExcitingStructure


def mock_species_files(tmpdir: py.path.local, species: List[str]):
    """
    Mock species files for test cases.

    :param tmpdir: Temporary directory, used by pytest
    :param species: List of species to mock
    """
    for element in species:
        file = tmpdir / element + ".xml"
        file.write("Arbitrary text")


@pytest.fixture
def xml_structure_H2He(tmpdir):
    """
    structure object initialised with a mock crystal, using mandatory arguments only.
    """
    cubic_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    arbitrary_atoms = [{'species': 'H', 'position': [0, 0, 0]},
                       {'species': 'H', 'position': [1, 0, 0]},
                       {'species': 'He', 'position': [2, 0, 0]}]
    mock_species_files(tmpdir, ["H", "He"])
    structure = ExcitingStructure(arbitrary_atoms, cubic_lattice, tmpdir)
    return structure.to_xml()


def test_class_exciting_structure_xml(tmpdir, xml_structure_H2He):
    """
    Test input XML attributes from an instance of ExcitingStructure.
    """
    assert xml_structure_H2He.tag == 'structure', 'XML root should be structure'
    assert xml_structure_H2He.keys() == ['speciespath'], 'structure defined to have only speciespath '
    assert xml_structure_H2He.get('speciespath') == str(tmpdir), 'species path set to tmpdir'


def test_class_exciting_structure_crystal_xml(tmpdir, xml_structure_H2He):
    """
    crystal subtree of structure.
    """
    elements = list(xml_structure_H2He)
    assert len(elements) == 3, 'Expect structure tree to have 3 sub-elements'

    crystal_xml = elements[0]
    assert crystal_xml.tag == "crystal", 'First subtree is crystal'
    assert crystal_xml.items() == [], "No required attributes in crystal."

    lattice_vectors = list(crystal_xml)
    assert len(lattice_vectors) == 3, "Always expect three lattice vectors"
    assert lattice_vectors[0].items() == [], 'Lattice vectors have no items'
    assert lattice_vectors[0].text == "1.0 0.0 0.0", "Lattice vector `a` differs from input"
    assert lattice_vectors[1].text == "0.0 1.0 0.0", "Lattice vector `b` differs from input"
    assert lattice_vectors[2].text == "0.0 0.0 1.0", "Lattice vector `c` differs from input"


def test_class_exciting_structure_species_xml(tmpdir, xml_structure_H2He):
    """
    species subtrees of structure.
    """
    elements = list(xml_structure_H2He)
    assert len(elements) == 3, 'Expect structure tree to have 3 sub-elements'

    species_h_xml = elements[1]
    assert species_h_xml.tag == "species", 'Second subtree is species'

    species_he_xml = elements[2]
    assert species_he_xml.tag == "species", 'Third subtree is species'

    assert species_h_xml.items() == [('speciesfile', 'H.xml')], 'species is inconsistent'
    assert species_he_xml.items() == [('speciesfile', 'He.xml')], 'species is inconsistent'

    atoms_h = list(species_h_xml)
    assert len(atoms_h) == 2, 'Number of H atoms is wrong'
    assert atoms_h[0].items() == [('coord', '0 0 0')], "Coordinate of first H differs to input"
    assert atoms_h[1].items() == [('coord', '1 0 0')], "Coordinate of second H differs to input"

    atoms_he = list(species_he_xml)
    assert len(atoms_he) == 1, 'Number of He atoms is wrong'
    assert atoms_he[0].items() == [('coord', '2 0 0')], "Coordinate of only He differs to input"


@pytest.fixture
def xml_structure_CdS(tmpdir):
    """
    structure object initialised with a mock crystal, using all atom properties.
    Optional atom attributes require the generic container, List[dict].
    """
    cubic_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    arbitrary_atoms = [{'species': 'Cd',
                        'position': [0.0, 0.0, 0.0],
                        'bfcmt': [1.0, 1.0, 1.0],
                        'lockxyz': [False, False, False],
                        'mommtfix': [2.0, 2.0, 2.0]},
                       {'species': 'S',
                        'position': [1.0, 0.0, 0.0]}
                       ]
    mock_species_files(tmpdir, ["Cd", "S"])
    structure = ExcitingStructure(arbitrary_atoms, cubic_lattice, str(tmpdir))
    return structure.to_xml()


def test_optional_atom_attributes_xml(tmpdir, xml_structure_CdS):
    """
    Test optional atom attributes are correctly set in XML tree.
    """
    assert xml_structure_CdS.tag == 'structure'
    assert xml_structure_CdS.keys() == ['speciespath'], 'structure defined to have only speciespath '
    assert xml_structure_CdS.get('speciespath') == str(tmpdir), 'speciespath should equal tmpdir'

    elements = list(xml_structure_CdS)
    assert len(elements) == 3, 'Expect structure tree to have 3 sub-elements'

    # Crystal
    crystal_xml = elements[0]
    assert crystal_xml.tag == "crystal", 'First subtree is crystal'
    assert crystal_xml.items() == [], "No required attributes in crystal."

    # Species
    species_cd_xml = elements[1]
    assert species_cd_xml.tag == "species", 'Second subtree is species'
    assert species_cd_xml.items() == [('speciesfile', 'Cd.xml')]

    species_s_xml = elements[2]
    assert species_s_xml.tag == "species", 'Third subtree is species'
    assert species_s_xml.items() == [('speciesfile', 'S.xml')]

    # Cd
    atoms_cd = list(species_cd_xml)
    assert len(atoms_cd) == 1, 'Wrong number of Cd atoms'
    assert set(atoms_cd[0].keys()) == {'coord', 'bfcmt', 'mommtfix', 'lockxyz'}, \
        'Cd contains all mandatory and optional atom properties'
    assert ('coord', '0.0 0.0 0.0') in atoms_cd[0].items()
    assert ('bfcmt', '1.0 1.0 1.0') in atoms_cd[0].items()
    assert ('mommtfix', '2.0 2.0 2.0') in atoms_cd[0].items()
    assert ('lockxyz', 'false false false') in atoms_cd[0].items()

    # S
    atoms_s = list(species_s_xml)
    assert len(atoms_s) == 1, 'Wrong number of S atoms'
    assert atoms_s[0].keys() == ['coord'], \
        'S only contains mandatory atom properties'
    assert atoms_s[0].items() == [('coord', '1.0 0.0 0.0')]


@pytest.fixture
def lattice_and_atoms_CdS(tmpdir):
    """
    structure object initialised with a mock crystal
    """
    cubic_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    arbitrary_atoms = [{'species': 'Cd', 'position': [0.0, 0.0, 0.0]},
                       {'species': 'S', 'position': [1.0, 0.0, 0.0]}]
    mock_species_files(tmpdir, ["Cd", "S"])
    return cubic_lattice, arbitrary_atoms


def test_optional_structure_attributes_xml(tmpdir, lattice_and_atoms_CdS):
    """
    Test optional structure attributes.
    """
    cubic_lattice, arbitrary_atoms = lattice_and_atoms_CdS
    structure_attributes = {'autormt': True,
                            'cartesian': False,
                            'epslat': 1.e-6,
                            'primcell': False,
                            'tshift': True
                            }
    structure = ExcitingStructure(arbitrary_atoms,
                                  cubic_lattice,
                                  str(tmpdir),
                                  structure_properties=structure_attributes)
    xml_structure = structure.to_xml()

    mandatory = ['speciespath']
    optional = list(structure_attributes)

    assert xml_structure.tag == 'structure'
    assert xml_structure.keys() == mandatory + optional, \
        'Should contain mandatory speciespath plus all optional attributes'
    assert xml_structure.get('speciespath') == str(tmpdir), 'species path should equal tmpdir'
    assert xml_structure.get('autormt') == 'true'
    assert xml_structure.get('cartesian') == 'false'
    assert xml_structure.get('epslat') == '1e-06'
    assert xml_structure.get('primcell') == 'false'
    assert xml_structure.get('tshift') == 'true'


def test_optional_crystal_attributes_xml(tmpdir, lattice_and_atoms_CdS):
    """
    Test optional crystal attributes.
    """
    cubic_lattice, arbitrary_atoms = lattice_and_atoms_CdS

    structure = ExcitingStructure(arbitrary_atoms,
                                  cubic_lattice,
                                  str(tmpdir),
                                  crystal_properties={'scale': 1.00, 'stretch': 1.00}
                                  )
    xml_structure = structure.to_xml()

    elements = list(xml_structure)
    assert len(elements) == 3, 'Number of sub-elements in structure tree'

    crystal_xml = elements[0]
    assert crystal_xml.tag == "crystal", 'First subtree is crystal'
    assert crystal_xml.keys() == ['scale', 'stretch'], 'Optional crystal properties'
    assert crystal_xml.get('scale') == '1.0', 'scale value inconsistent with input'
    assert crystal_xml.get('stretch') == '1.0', 'stretch value inconsistent with input'


def test_optional_species_attributes_xml(tmpdir, lattice_and_atoms_CdS):
    """
    Test optional species attributes.
    """
    cubic_lattice, arbitrary_atoms = lattice_and_atoms_CdS
    species_attributes = {'Cd': {'rmt': 3.0}, 'S': {'rmt': 4.0}}

    structure = ExcitingStructure(arbitrary_atoms,
                                  cubic_lattice,
                                  str(tmpdir),
                                  species_properties=species_attributes)
    xml_structure = structure.to_xml()

    elements = list(xml_structure)
    assert len(elements) == 3, 'Number of sub-elements in structure tree'

    species_cd_xml = elements[1]
    assert species_cd_xml.tag == "species", 'Second subtree is species'

    species_s_xml = elements[2]
    assert species_s_xml.tag == "species", 'Third subtree is species'

    assert species_cd_xml.keys() == ['speciesfile', 'rmt'], "species attributes differ from expected"
    assert species_cd_xml.get('speciesfile') == 'Cd.xml', 'speciesfile differs from expected'
    assert species_cd_xml.get('rmt') == '3.0', 'Cd muffin tin radius differs from input'

    assert species_s_xml.keys() == ['speciesfile', 'rmt'], "species attributes differ from expected"
    assert species_s_xml.get('speciesfile') == 'S.xml', 'speciesfile differs from expected'
    assert species_s_xml.get('rmt') == '4.0', 'S muffin tin radius differs from input'


@pytest.fixture
def lattice_and_atoms_H20(tmpdir):
    """
    structure object initialised with a H20 molecule in a big box (angstrom)
    """
    cubic_lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
    arbitrary_atoms = [{'species': 'H', 'position': [0.00000, 0.75545, -0.47116]},
                       {'species': 'O', 'position': [0.00000, 0.00000, 0.11779]},
                       {'species': 'H', 'position': [0.00000, 0.75545, -0.47116]}]
    mock_species_files(tmpdir, ["H", "O"])
    return cubic_lattice, arbitrary_atoms


def test_group_atoms_by_species(tmpdir, lattice_and_atoms_H20):
    """
    Test group_atoms_by_species method.
    """
    cubic_lattice, atoms = lattice_and_atoms_H20
    structure = ExcitingStructure(atoms, cubic_lattice, tmpdir)
    assert structure.species == ['H', 'O', 'H'], 'Species list differs from lattice_and_atoms_H20'

    indices = structure._group_atoms_by_species()
    assert set(indices.keys()) == {'H', 'O'}, 'List unique species in structure'
    assert indices['H'] == [0, 2], "Expect atoms 0 and 2 to be H"
    assert indices['O'] == [1], "Expect atom 1 to be O"

    hydrogen = [structure.species[i] for i in indices['H']]
    oxygen = [structure.species[i] for i in indices['O']]
    assert hydrogen == ["H", "H"], "Expect list to only contain H symbols"
    assert oxygen == ["O"], "Expect list to contain only one O symbol"


def test_class_exciting_structure_ase(tmpdir, lattice_and_atoms_H20):
    """
    Test the ASE Atoms object gets used correctly by the ExcitingStructure constructor.
    """
    if "ase" not in sys.modules:
        # ASE not available, so do not run
        return

    cubic_cell = np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]])
    pos = [[0.00000, 0.75545, -0.47116],
           [0.00000, 0.00000, 0.11779],
           [0.00000, 0.75545, -0.47116]]
    ase_atoms = ase.atoms.Atoms(symbols=["H", "O", "H"], positions=pos, cell=cubic_cell)
    mock_species_files(tmpdir, ["H", "O", "H"])

    structure = ExcitingStructure(ase_atoms, species_path=tmpdir)

    assert structure.species == ["H", "O", "H"]
    assert np.allclose(structure.lattice,
                       [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]), \
        'Expect lattice vectors to match input values'

    assert np.allclose(structure.positions, pos), 'Expect positions to match input values.'


def test_check_for_species_files_not_exist_dir():
    """
    Ensure an error raised when species files directory does not exist.
    """
    cubic_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    arbitrary_atoms = [{'species': 'Cd', 'position': [0.0, 0.0, 0.0]},
                       {'species': 'S', 'position': [1.0, 0.0, 0.0]}]

    with pytest.raises(StopIteration) as error:
        ExcitingStructure(arbitrary_atoms, cubic_lattice, species_path="./not_existing_dir")

    assert error.value.args[0] == \
           "species_path does not exist or is empty: ./not_existing_dir"


def test_check_for_species_files_no_species_files(tmpdir):
    """
    Test that an error is raised when species files don't exist.
    Achieved by not mocking the species files in this test.
    """
    cubic_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    arbitrary_atoms = [{'species': 'Cd', 'position': [0.0, 0.0, 0.0]},
                       {'species': 'S', 'position': [1.0, 0.0, 0.0]}]

    with pytest.raises(FileNotFoundError) as error:
        ExcitingStructure(arbitrary_atoms, cubic_lattice, species_path=tmpdir)

    assert error.value.args[0] == f"species_path directory does not contain species files: {tmpdir}"


def test_check_for_species_files_missing_species(tmpdir):
    """
    Ensure error raised when a required species file is missing.
    """
    cubic_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    arbitrary_atoms = [{'species': 'H', 'position': [0.0, 0.0, 0.0]}]
    mock_species_files(tmpdir, ['He'])

    with pytest.raises(FileNotFoundError) as error:
        ExcitingStructure(arbitrary_atoms, cubic_lattice, species_path=tmpdir)

    assert error.value.args[0] == f"Species files, ['H.xml'], are missing from: {tmpdir}"


def test_initialise_with_bad_chemical_symbol(tmpdir):
    """
    Ensure error raised when initialising with a species that is not a
    chemical element.
    """
    cubic_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    arbitrary_atoms = [{'species': 'Xx', 'position': [1, 0, 0]}]

    with pytest.raises(ValueError) as error:
        ExcitingStructure(arbitrary_atoms, cubic_lattice, tmpdir)

    assert error.value.args[0] == "Species input keys are not valid: {'xx'}"
