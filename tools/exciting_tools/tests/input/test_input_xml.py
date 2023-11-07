"""Test composition of an exciting input XML.

TODO(Fab/Alex/Dan) Issue 117. Would be nice to assert that the output is valid
    XML * https://lxml.de/validation.html
Also see: https://xmlschema.readthedocs.io/en/latest/usage.html#xsd-declarations

NOTE:
All attribute tests should assert on the XML tree content's as the attribute
order is not preserved by the ElementTree.tostring method. Elements appear to
be fine.
"""

import pytest

from excitingtools.input.input_classes import ExcitingGroundStateInput, ExcitingXSInput, ExcitingInputXML
from excitingtools.input.structure import ExcitingStructure


@pytest.fixture
def exciting_structure() -> ExcitingStructure:
    """ Initialise an exciting structure. """
    cubic_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    arbitrary_atoms = [{'species': 'Li', 'position': [0.0, 0.0, 0.0]},
                       {'species': 'Li', 'position': [1.0, 0.0, 0.0]},
                       {'species': 'F', 'position': [2.0, 0.0, 0.0]}]

    return ExcitingStructure(arbitrary_atoms, cubic_lattice, '.')


@pytest.fixture
def exciting_input_xml(exciting_structure: ExcitingStructure) -> ExcitingInputXML:
    """ Initialises a complete input file. """
    ground_state = ExcitingGroundStateInput(
        rgkmax=8.0,
        do="fromscratch",
        ngridk=[6, 6, 6],
        xctype="GGA_PBE_SOL",
        vkloff=[0, 0, 0],
        tforce=True,
        nosource=False
    )

    bse_attributes = {'bsetype': 'singlet', 'xas': True}
    energywindow_attributes = {'intv': [5.8, 8.3], 'points': 5000}
    screening_attributes = {'screentype': 'full', 'nempty': 15}
    plan_input = ['screen', 'bse']
    qpointset_input = [[0, 0, 0], [0.5, 0.5, 0.5]]
    xs = ExcitingXSInput(xstype="BSE", broad=0.32, ngridk=[8, 8, 8],
                         BSE=bse_attributes,
                         energywindow=energywindow_attributes,
                         screening=screening_attributes,
                         qpointset=qpointset_input,
                         plan=plan_input)

    return ExcitingInputXML(sharedfs=True, structure=exciting_structure, title='Test Case', groundstate=ground_state,
                            xs=xs)


def test_exciting_input_xml_structure_and_gs_and_xs(exciting_input_xml: ExcitingInputXML):
    """Test the XML created for a ground state input is valid.
    Test SubTree composition using only mandatory attributes for each XML subtree.
    """
    input_xml_tree = exciting_input_xml.to_xml()

    assert input_xml_tree.tag == 'input'
    assert input_xml_tree.keys() == ['sharedfs']

    subelements = list(input_xml_tree)
    assert len(subelements) == 4

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

    xs_xml = subelements[3]
    assert xs_xml.tag == 'xs'
    assert set(xs_xml.keys()) == {'broad', 'ngridk', 'xstype'}
    assert xs_xml.get('broad') == '0.32'
    assert xs_xml.get('ngridk') == '8 8 8'
    assert xs_xml.get('xstype') == 'BSE'

    xs_subelements = list(xs_xml)
    assert len(xs_subelements) == 5
    valid_tags = {"screening", "BSE", "energywindow", "qpointset", "plan"}
    assert valid_tags == set(xs_subelement.tag for xs_subelement in xs_subelements)

    screening_xml = xs_xml.find("screening")
    assert screening_xml.tag == "screening"
    assert screening_xml.keys() == ['screentype', 'nempty']
    assert screening_xml.get('screentype') == 'full'
    assert screening_xml.get('nempty') == '15'

    bse_xml = xs_xml.find("BSE")
    assert bse_xml.tag == 'BSE'
    assert bse_xml.keys() == ['bsetype', 'xas']
    assert bse_xml.get('bsetype') == 'singlet'
    assert bse_xml.get('xas') == 'true'

    energywindow_xml = xs_xml.find("energywindow")
    assert energywindow_xml.tag == "energywindow"
    assert energywindow_xml.keys() == ['intv', 'points']
    assert energywindow_xml.get('intv') == '5.8 8.3'
    assert energywindow_xml.get('points') == '5000'

    qpointset_xml = xs_xml.find("qpointset")
    assert qpointset_xml.tag == "qpointset"
    assert qpointset_xml.items() == []
    qpoints = list(qpointset_xml)
    assert len(qpoints) == 2
    assert qpoints[0].tag == 'qpoint'
    assert qpoints[0].items() == []
    valid_qpoints = {'0 0 0', '0.5 0.5 0.5'}
    assert qpoints[0].text in valid_qpoints
    valid_qpoints.discard(qpoints[0].text)
    assert qpoints[1].text in valid_qpoints

    plan_xml = xs_xml.find("plan")
    assert plan_xml.tag == "plan"
    assert plan_xml.items() == []
    doonlys = list(plan_xml)
    assert len(doonlys) == 2
    assert doonlys[0].tag == 'doonly'
    assert doonlys[0].items() == [('task', 'screen')]
    assert doonlys[1].tag == 'doonly'
    assert doonlys[1].items() == [('task', 'bse')]


def test_attribute_modification(exciting_input_xml: ExcitingInputXML):
    """Test the XML created for a ground state input is valid.
    Test SubTree composition using only mandatory attributes for each XML subtree.
    """
    exciting_input_xml.set_title("New Test Case")
    exciting_input_xml.structure.crystal_properties.scale = 2.3
    exciting_input_xml.groundstate.rgkmax = 9.0
    exciting_input_xml.xs.energywindow.points = 4000
    input_xml_tree = exciting_input_xml.to_xml()

    subelements = list(input_xml_tree)
    assert len(subelements) == 4

    title_xml = subelements[0]
    assert title_xml.tag == 'title'
    assert title_xml.text == 'New Test Case'

    structure_xml = subelements[1]
    assert structure_xml[0].get("scale") == "2.3"

    groundstate_xml = subelements[2]
    assert groundstate_xml.get('rgkmax') == "9.0"

    xs_xml = subelements[3]
    xs_subelements = list(xs_xml)
    assert len(xs_subelements) == 5

    energywindow_xml = xs_xml.find("energywindow")
    assert energywindow_xml.get('points') == '4000'


def test_as_dict(exciting_input_xml: ExcitingInputXML, mock_env_jobflow_missing):
    dict_representation = exciting_input_xml.as_dict()
    assert set(dict_representation.keys()) == {"xml_string"}
    # check only that the xml string starts with the correct first lines:
    assert dict_representation["xml_string"].startswith('<?xml version="1.0" ?>\n<input sharedfs="true">\n\t'
                                                        '<title>Test Case</title>\n\t<structure')


def test_as_dict_jobflow(exciting_input_xml: ExcitingInputXML, mock_env_jobflow):
    dict_representation = exciting_input_xml.as_dict()
    xml_string = dict_representation.pop("xml_string")
    assert dict_representation == {'@class': 'ExcitingInputXML',
                                   '@module': 'excitingtools.input.input_classes'}
    # check only that the xml string starts with the correct first lines:
    assert xml_string.startswith('<?xml version="1.0" ?>\n<input sharedfs="true">\n\t'
                                 '<title>Test Case</title>\n\t<structure')


def test_from_dict(exciting_input_xml: ExcitingInputXML, mock_env_jobflow_missing):
    new_input_xml = ExcitingInputXML.from_dict(exciting_input_xml.as_dict())
    assert new_input_xml.title.title == "Test Case"  # pylint: disable=no-member
    assert new_input_xml.groundstate.ngridk == [6, 6, 6]  # pylint: disable=no-member


def test_from_dict_jobflow(exciting_input_xml: ExcitingInputXML, mock_env_jobflow):
    new_input_xml = ExcitingInputXML.from_dict(exciting_input_xml.as_dict())
    assert new_input_xml.title.title == "Test Case"  # pylint: disable=no-member


def test_missing_structure():
    with pytest.raises(ValueError, match="Missing mandatory arguments: {'structure'}"):
        ExcitingInputXML(title="Test Case", groundstate=ExcitingGroundStateInput())


def test_from_gs_dict(exciting_structure):
    input_xml = ExcitingInputXML(structure=exciting_structure, title="Test Case",
                                 groundstate={"rgkmax": 7.0})
    assert input_xml.title.title == "Test Case"  # pylint: disable=no-member
    assert input_xml.groundstate.rgkmax == 7.0  # pylint: disable=no-member
    assert input_xml.structure.speciespath == "."  # pylint: disable=no-member
