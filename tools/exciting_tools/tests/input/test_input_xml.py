"""Test composition of an exciting input XML.

TODO(Fab/Alex/Dan) Issue 117. Would be nice to assert that the output is valid
    XML * https://lxml.de/validation.html
"""
from typing import List

import py

from excitingtools.input.input_xml import exciting_input_xml
from excitingtools.input.structure import ExcitingStructure
from excitingtools.input.ground_state import ExcitingGroundStateInput


def mock_species_files(species_dir: py.path.local, species: List[str]):
    """
    Mock species files for test cases.

    :param species_dir: Temporary directory, used by pytest
    :param species: List of species to mock
    """
    for x in species:
        file = species_dir / x + ".xml"
        file.write("Arbitrary text")


def test_exciting_input_xml_structure_and_gs(tmpdir):
    """Test the XML created for a ground state input is valid.
    Test SubTree composition using only mandatory attributes for each XML subtree.
    """
    # Structure
    cubic_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    arbitrary_atoms = [{'species': 'Li', 'position': [0.0, 0.0, 0.0]},
                       {'species': 'Li', 'position': [1.0, 0.0, 0.0]},
                       {'species': 'F', 'position': [2.0, 0.0, 0.0]}]

    mock_species_files(tmpdir, ["Li", "F"])
    structure = ExcitingStructure(arbitrary_atoms, cubic_lattice, tmpdir)

    ground_state = ExcitingGroundStateInput(
        rgkmax=8.0,
        do="fromscratch",
        ngridk=[6, 6, 6],
        xctype="GGA_PBE_SOL",
        vkloff=[0, 0, 0],
        tforce=True,
        nosource=False
        )

    input_xml_tree = exciting_input_xml(structure, title='Test Case', groundstate=ground_state)

    assert input_xml_tree.tag == 'input'
    assert input_xml_tree.keys() == []

    subelements = list(input_xml_tree)
    assert len(subelements) == 3

    title_xml = subelements[0]
    assert title_xml.tag == 'title'
    assert title_xml.keys() == []
    assert title_xml.text == 'Test Case'

    structure_xml = subelements[1]
    assert structure_xml.tag == 'structure'
    assert structure_xml.keys() == ['speciespath']
    assert len(list(structure_xml)) == 3

    groundstate_xml = subelements[2]
    assert groundstate_xml.tag == 'groundstate'
    assert groundstate_xml.text == ' '
    assert groundstate_xml.keys() == \
           ['rgkmax', 'do', 'ngridk', 'xctype', 'vkloff', 'tforce', 'nosource']
    assert groundstate_xml.get('rgkmax') == "8.0"
    assert groundstate_xml.get('do') == "fromscratch"
    assert groundstate_xml.get('ngridk') == "6 6 6"
    assert groundstate_xml.get('xctype') == "GGA_PBE_SOL"
    assert groundstate_xml.get('vkloff') == "0 0 0"
    assert groundstate_xml.get('tforce') == "true"
    assert groundstate_xml.get('nosource') == "false"
