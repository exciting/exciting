"""Test ExcitingStructure, python API that generates exciting's structure XML.

NOTE:
All attribute tests should assert on the XML tree content, as the attribute
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

import numpy as np
import pytest

from excitingtools.input.bandstructure import get_bandstructure_input_from_exciting_structure
from excitingtools.input.structure import ExcitingStructure


@pytest.fixture
def xml_structure_H2He():
    """
    structure object initialised with a mock crystal, using mandatory arguments only.
    """
    cubic_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    arbitrary_atoms = [{'species': 'H', 'position': [0, 0, 0]},
                       {'species': 'H', 'position': [1, 0, 0]},
                       {'species': 'He', 'position': [2, 0, 0]}]
    structure = ExcitingStructure(arbitrary_atoms, cubic_lattice, './')
    return structure.to_xml()


def test_class_exciting_structure_xml(xml_structure_H2He):
    """
    Test input XML attributes from an instance of ExcitingStructure.
    """
    assert xml_structure_H2He.tag == 'structure', 'XML root should be structure'
    assert xml_structure_H2He.keys() == ['speciespath'], 'structure defined to have only speciespath '
    assert xml_structure_H2He.get('speciespath') == './', 'species path set to ./'


def test_class_exciting_structure_crystal_xml(xml_structure_H2He):
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


def test_class_exciting_structure_species_xml(xml_structure_H2He):
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
def xml_structure_CdS():
    """
    structure object initialised with a mock crystal, using all atom properties.
    Optional atom attributes require the generic container, List[dict].
    """
    cubic_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    arbitrary_atoms = [{
        'species': 'Cd',
        'position': [0.0, 0.0, 0.0],
        'bfcmt': [1.0, 1.0, 1.0],
        'lockxyz': [False, False, False],
        'mommtfix': [2.0, 2.0, 2.0]
    }, {'species': 'S', 'position': [1.0, 0.0, 0.0]}]
    structure = ExcitingStructure(arbitrary_atoms, cubic_lattice, './')
    return structure.to_xml()


def test_optional_atom_attributes_xml(xml_structure_CdS):
    """
    Test optional atom attributes are correctly set in XML tree.
    """
    assert xml_structure_CdS.tag == 'structure'
    assert xml_structure_CdS.keys() == ['speciespath'], 'structure defined to have only speciespath '
    assert xml_structure_CdS.get('speciespath') == './', 'speciespath set to be ./'

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
def lattice_and_atoms_CdS():
    """
    structure object initialised with a mock crystal
    """
    cubic_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    arbitrary_atoms = [{'species': 'Cd', 'position': [0.0, 0.0, 0.0]},
                       {'species': 'S', 'position': [1.0, 0.0, 0.0]}]
    return cubic_lattice, arbitrary_atoms


def test_optional_structure_attributes_xml(lattice_and_atoms_CdS):
    """
    Test optional structure attributes.
    """
    cubic_lattice, arbitrary_atoms = lattice_and_atoms_CdS
    structure_attributes = {
        'autormt': True, 'cartesian': False, 'epslat': 1.e-6, 'primcell': False, 'tshift': True
    }
    structure = ExcitingStructure(arbitrary_atoms, cubic_lattice, './', **structure_attributes)
    xml_structure = structure.to_xml()

    mandatory = {'speciespath'}
    optional = set(structure_attributes)

    assert xml_structure.tag == 'structure'
    assert set(xml_structure.keys()) == mandatory | optional, \
        'Should contain mandatory speciespath plus all optional attributes'
    assert xml_structure.get('speciespath') == './', 'species path should be ./'
    assert xml_structure.get('autormt') == 'true'
    assert xml_structure.get('cartesian') == 'false'
    assert xml_structure.get('epslat') == '1e-06'
    assert xml_structure.get('primcell') == 'false'
    assert xml_structure.get('tshift') == 'true'


def test_optional_crystal_attributes_xml(lattice_and_atoms_CdS):
    """
    Test optional crystal attributes.
    """
    cubic_lattice, arbitrary_atoms = lattice_and_atoms_CdS

    structure = ExcitingStructure(
        arbitrary_atoms,
        cubic_lattice,
        './',
        crystal_properties={'scale': 1.00, 'stretch': [1.00, 1.00, 1.00]}
    )
    xml_structure = structure.to_xml()

    elements = list(xml_structure)
    assert len(elements) == 3, 'Number of sub-elements in structure tree'

    crystal_xml = elements[0]
    assert crystal_xml.tag == "crystal", 'First subtree is crystal'
    assert crystal_xml.keys() == ['scale', 'stretch'], 'Optional crystal properties'
    assert crystal_xml.get('scale') == '1.0', 'scale value inconsistent with input'
    assert crystal_xml.get('stretch') == '1.0 1.0 1.0', 'stretch value inconsistent with input'


def test_optional_species_attributes_xml(lattice_and_atoms_CdS):
    """
    Test optional species attributes.
    """
    cubic_lattice, arbitrary_atoms = lattice_and_atoms_CdS
    species_attributes = {'Cd': {'rmt': 3.0, "LDAplusU": {"J": 1.5, "U": 2.4, "l": 2}},
                          'S': {'rmt': 4.0, "dfthalfparam": {"ampl": 1.2, "cut": 1.9, "exponent": 5,
                                                             "shell": [{"ionization": 0.8, "number": 1}]}}}

    structure = ExcitingStructure(
        arbitrary_atoms, cubic_lattice, './', species_properties=species_attributes
    )
    xml_structure = structure.to_xml()

    elements = list(xml_structure)
    assert len(elements) == 3, 'Number of sub-elements in structure tree'

    species_cd_xml = elements[1]
    assert species_cd_xml.tag == "species", 'Second subtree is species'

    species_s_xml = elements[2]
    assert species_s_xml.tag == "species", 'Third subtree is species'

    assert set(species_cd_xml.keys()) == {'speciesfile', 'rmt'}, "species attributes differ from expected"
    assert species_cd_xml.get('speciesfile') == 'Cd.xml', 'speciesfile differs from expected'
    assert species_cd_xml.get('rmt') == '3.0', 'Cd muffin tin radius differs from input'

    species_cd_elements = list(species_cd_xml)
    assert len(species_cd_elements) == 2
    ldaplusu_xml = species_cd_elements[0]
    assert ldaplusu_xml.tag == "LDAplusU"

    assert set(ldaplusu_xml.keys()) == {'J', 'U', 'l'}
    assert ldaplusu_xml.get("J") == "1.5"
    assert ldaplusu_xml.get("U") == "2.4"
    assert ldaplusu_xml.get("l") == "2"

    assert set(species_s_xml.keys()) == {'speciesfile', 'rmt'}, "species attributes differ from expected"
    assert species_s_xml.get('speciesfile') == 'S.xml', 'speciesfile differs from expected'
    assert species_s_xml.get('rmt') == '4.0', 'S muffin tin radius differs from input'

    species_s_elements = list(species_s_xml)
    assert len(species_s_elements) == 2
    dfthalfparam_xml = species_s_elements[0]
    assert dfthalfparam_xml.tag == "dfthalfparam"

    assert set(dfthalfparam_xml.keys()) == {'ampl', 'cut', 'exponent'}
    assert dfthalfparam_xml.get("ampl") == "1.2"
    assert dfthalfparam_xml.get("cut") == "1.9"
    assert dfthalfparam_xml.get("exponent") == "5"

    dfthalfparam_elements = list(dfthalfparam_xml)
    assert len(dfthalfparam_elements) == 1
    shell_xml = dfthalfparam_elements[0]
    assert shell_xml.tag == "shell"

    assert set(shell_xml.keys()) == {'ionization', 'number'}
    assert shell_xml.get('ionization') == '0.8'
    assert shell_xml.get('number') == '1'


ref_dict = {'xml_string': '<structure speciespath="./"> <crystal> <basevect>1.0 0.0 0.0</basevect>'
                          '<basevect>0.0 1.0 0.0</basevect><basevect>0.0 0.0 1.0</basevect></crystal>'
                          '<species speciesfile="Cd.xml"> <atom coord="0.0 0.0 0.0"> </atom></species>'
                          '<species speciesfile="S.xml"> <atom coord="1.0 0.0 0.0"> </atom></species></structure>'}


def test_as_dict(lattice_and_atoms_CdS, mock_env_jobflow_missing):
    cubic_lattice, arbitrary_atoms = lattice_and_atoms_CdS

    structure = ExcitingStructure(arbitrary_atoms, cubic_lattice, './')

    assert structure.as_dict() == ref_dict, 'expected different dict representation'


def test_as_dict_jobflow(lattice_and_atoms_CdS, mock_env_jobflow):
    cubic_lattice, arbitrary_atoms = lattice_and_atoms_CdS

    structure = ExcitingStructure(arbitrary_atoms, cubic_lattice, './')

    assert structure.as_dict() == {**ref_dict, '@class': 'ExcitingStructure',
                                   '@module': 'excitingtools.input.structure'}, \
        'expected different dict representation'


def test_from_dict(lattice_and_atoms_CdS):
    cubic_lattice, arbitrary_atoms = lattice_and_atoms_CdS
    structure = ExcitingStructure.from_dict(ref_dict)

    assert np.allclose(structure.lattice, np.array(cubic_lattice))
    assert structure.species == [d['species'] for d in arbitrary_atoms]
    assert structure.positions == [d['position'] for d in arbitrary_atoms]
    assert structure.speciespath == './'  # pylint: disable=no-member


@pytest.fixture
def lattice_and_atoms_H20():
    """
    H20 molecule in a big box (angstrom)
    """
    cubic_lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
    atoms = [{'species': 'H', 'position': [0.00000, 0.75545, -0.47116]},
             {'species': 'O', 'position': [0.00000, 0.00000, 0.11779]},
             {'species': 'H', 'position': [0.00000, 0.75545, -0.47116]}]
    return cubic_lattice, atoms


def test_group_atoms_by_species(lattice_and_atoms_H20):
    """
    Test group_atoms_by_species method.
    """
    cubic_lattice, atoms = lattice_and_atoms_H20
    structure = ExcitingStructure(atoms, cubic_lattice, './')
    assert structure.species == ['H', 'O', 'H'], 'Species list differs from lattice_and_atoms_H20'

    indices = dict(structure._group_atoms_by_species())
    assert set(indices.keys()) == {'H', 'O'}, 'List unique species in structure'
    assert indices['H'] == [0, 2], "Expect atoms 0 and 2 to be H"
    assert indices['O'] == [1], "Expect atom 1 to be O"

    hydrogen = [structure.species[i] for i in indices['H']]
    oxygen = [structure.species[i] for i in indices['O']]
    assert hydrogen == ["H", "H"], "Expect list to only contain H symbols"
    assert oxygen == ["O"], "Expect list to contain only one O symbol"


@pytest.fixture
def ase_atoms_H20(lattice_and_atoms_H20):
    """
    H20 molecule in a big box (angstrom), in ASE Atoms()
    Converts a List[dict] to ase.atoms.Atoms.
    """
    ase = pytest.importorskip("ase")
    lattice, atoms = lattice_and_atoms_H20
    symbols = [atom['species'] for atom in atoms]
    cubic_cell = np.asarray(lattice)
    positions = [atom['position'] for atom in atoms]
    return ase.Atoms(symbols=symbols, positions=positions, cell=cubic_cell, pbc=True)


def test_class_exciting_structure_ase(ase_atoms_H20):
    """
    Test the ASE Atoms object gets used correctly by the ExcitingStructure constructor.
    """
    atoms = ase_atoms_H20
    structure = ExcitingStructure(atoms, species_path='./')

    assert structure.species == ["H", "O", "H"]
    assert np.allclose(structure.lattice,
                       [[18.897261246257703, 0.0, 0.0],
                        [0.0, 18.897261246257703, 0.0],
                        [0.0, 0.0, 18.897261246257703]]), \
        'Expect lattice vectors to match input values'

    assert np.allclose(structure.positions, atoms.get_scaled_positions()), 'Expect positions to match input values.'

    # TODO(Alex) Issue 117. Compare xml_structure built with and without ASE - should be consistent
    # This just confirms the XML tree is built, not that it is correct.
    xml_structure = structure.to_xml()
    assert list(xml_structure.keys()) == ['speciespath'], 'Only expect speciespath in structure xml keys'


def test_using_non_standard_species_symbol():
    cubic_lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
    atoms = [{'species': 'C_molecule', 'position': [0.00000, 0.75545, -0.47116]}]
    structure = ExcitingStructure(atoms, cubic_lattice)

    structure_xml = structure.to_xml()
    assert structure_xml.tag == 'structure', 'XML root should be structure'
    assert structure_xml.keys() == ['speciespath'], 'structure defined to have only speciespath '
    assert structure_xml.get('speciespath') == './', 'species path set to ./'

    elements = list(structure_xml)
    assert len(elements) == 2, 'Expect structure tree to have 2 sub-elements'

    species_c_xml = elements[1]
    assert species_c_xml.tag == "species", 'Second subtree is species'

    assert species_c_xml.items() == [('speciesfile', 'C_molecule.xml')], 'species is inconsistent'

    atoms_h = list(species_c_xml)
    assert len(atoms_h) == 1
    assert atoms_h[0].items() == [('coord', '0.0 0.75545 -0.47116')]


def test_get_full_lattice(lattice_and_atoms_CdS):
    cubic_lattice, arbitrary_atoms = lattice_and_atoms_CdS
    structure = ExcitingStructure(arbitrary_atoms, cubic_lattice, './',
                                  crystal_properties={'scale': 1.50, 'stretch': [2.00, 1.00, 3.00]})
    ref_lattice = np.array([[3, 0, 0], [0, 1.5, 0], [0, 0, 4.5]])
    assert np.allclose(structure.get_lattice(), ref_lattice)


def test_get_bandstructure_input_from_exciting_structure(lattice_and_atoms_H20):
    pytest.importorskip("ase")
    cubic_lattice, atoms = lattice_and_atoms_H20
    structure = ExcitingStructure(atoms, cubic_lattice, './')
    bandstructure = get_bandstructure_input_from_exciting_structure(structure)
    bs_xml = bandstructure.to_xml()

    assert bs_xml.tag == "bandstructure", 'Root tag should be "bandstructure"'

    plot1d_xml = bs_xml.find("plot1d")
    assert plot1d_xml is not None, 'Missing "plot1d" subtree'

    path_xml = plot1d_xml.find("path")
    assert path_xml is not None, 'Missing "path" subtree'
    assert path_xml.get("steps") == "100", 'Invalid value for "steps" attribute'

    assert len(list(path_xml)) == 8
    point1 = path_xml[0]
    assert point1.get("coord") == "0.0 0.0 0.0", 'Invalid value for "coord" attribute of point 1'
    assert point1.get("label") == "G", 'Invalid value for "label" attribute of point 1'

    point2 = path_xml[1]
    assert point2.get("coord") == "0.0 0.5 0.0", 'Invalid value for "coord" attribute of point 2'
    assert point2.get("label") == "X", 'Invalid value for "label" attribute of point 2'

    point3 = path_xml[2]
    assert point3.get("coord") == "0.5 0.5 0.0", 'Invalid value for "coord" attribute of point 3'
    assert point3.get("label") == "M", 'Invalid value for "label" attribute of point 3'

    point4 = path_xml[3]
    assert point4.get("coord") == "0.0 0.0 0.0", 'Invalid value for "coord" attribute of point 4'
    assert point4.get("label") == "G", 'Invalid value for "label" attribute of point 4'

    point5 = path_xml[4]
    assert point5.get("coord") == "0.5 0.5 0.5", 'Invalid value for "coord" attribute of point 5'
    assert point5.get("label") == "R", 'Invalid value for "label" attribute of point 5'

    point6 = path_xml[5]
    assert point6.get("coord") == "0.0 0.5 0.0", 'Invalid value for "coord" attribute of point 6'
    assert point6.get("label") == "X", 'Invalid value for "label" attribute of point 6'
    assert point6.get("breakafter") == "true", 'Invalid value for "breakafter" attribute of point 6'

    point7 = path_xml[6]
    assert point7.get("coord") == "0.5 0.5 0.0", 'Invalid value for "coord" attribute of point 7'
    assert point7.get("label") == "M", 'Invalid value for "label" attribute of point 7'

    point8 = path_xml[7]
    assert point8.get("coord") == "0.5 0.5 0.5", 'Invalid value for "coord" attribute of point 8'
    assert point8.get("label") == "R", 'Invalid value for "label" attribute of point 8'
