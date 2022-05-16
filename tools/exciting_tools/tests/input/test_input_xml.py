"""Test composition of an exciting input XML.

TODO(Fab/Alex/Dan) Issue 117. Would be nice to assert that the output is valid
    XML * https://lxml.de/validation.html
Also see: https://xmlschema.readthedocs.io/en/latest/usage.html#xsd-declarations
"""
from excitingtools.input.input_xml import exciting_input_xml
from excitingtools.input.structure import ExcitingStructure
from excitingtools.input.ground_state import ExcitingGroundStateInput
from excitingtools.input.xs import ExcitingXSInput


def test_exciting_input_xml_structure_and_gs_and_xs():
    """Test the XML created for a ground state input is valid.
    Test SubTree composition using only mandatory attributes for each XML subtree.
    """
    # Structure
    cubic_lattice = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    arbitrary_atoms = [{'species': 'Li', 'position': [0.0, 0.0, 0.0]},
                       {'species': 'Li', 'position': [1.0, 0.0, 0.0]},
                       {'species': 'F', 'position': [2.0, 0.0, 0.0]}]

    structure = ExcitingStructure(arbitrary_atoms, cubic_lattice, '.')

    ground_state = ExcitingGroundStateInput(
        rgkmax=8.0,
        do="fromscratch",
        ngridk=[6, 6, 6],
        xctype="GGA_PBE_SOL",
        vkloff=[0, 0, 0],
        tforce=True,
        nosource=False
        )

    xs_attributes = {'broad': 0.32, 'ngridk': [8, 8, 8]}
    bse_attributes = {'bsetype': 'singlet', 'xas': True}
    energywindow_attributes = {'intv': [5.8, 8.3], 'points': 5000}
    screening_attributes = {'screentype': 'full', 'nempty': 15}
    plan_input = ['screen', 'bse']
    qpointset_input = [[0, 0, 0], [0.5, 0.5, 0.5]]
    xs = ExcitingXSInput("BSE", xs=xs_attributes,
                         BSE=bse_attributes,
                         energywindow=energywindow_attributes,
                         screening=screening_attributes,
                         qpointset=qpointset_input,
                         plan=plan_input)

    input_xml_tree = exciting_input_xml(
        structure, title='Test Case', groundstate=ground_state, xs=xs)

    assert input_xml_tree.tag == 'input'
    assert input_xml_tree.keys() == []

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
    try:
        assert xs_xml.keys() == ['broad', 'ngridk', 'xstype']
    except AssertionError:
        assert xs_xml.keys() == ['xstype', 'broad', 'ngridk']
    assert xs_xml.get('broad') == '0.32'
    assert xs_xml.get('ngridk') == '8 8 8'
    assert xs_xml.get('xstype') == 'BSE'

    xs_subelements = list(xs_xml)
    assert len(xs_subelements) == 5

    screening_xml = xs_subelements[0]
    assert screening_xml.tag == "screening"
    assert screening_xml.keys() == ['screentype', 'nempty']
    assert screening_xml.get('screentype') == 'full'
    assert screening_xml.get('nempty') == '15'

    bse_xml = xs_subelements[1]
    assert bse_xml.tag == 'BSE'
    assert bse_xml.keys() == ['bsetype', 'xas']
    assert bse_xml.get('bsetype') == 'singlet'
    assert bse_xml.get('xas') == 'true'

    energywindow_xml = xs_subelements[2]
    assert energywindow_xml.tag == "energywindow"
    assert energywindow_xml.keys() == ['intv', 'points']
    assert energywindow_xml.get('intv') == '5.8 8.3'
    assert energywindow_xml.get('points') == '5000'

    qpointset_xml = xs_subelements[3]
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

    plan_xml = xs_subelements[4]
    assert plan_xml.tag == "plan"
    assert plan_xml.items() == []
    doonlys = list(plan_xml)
    assert len(doonlys) == 2
    assert doonlys[0].tag == 'doonly'
    assert doonlys[0].items() == [('task', 'screen')]
    assert doonlys[1].tag == 'doonly'
    assert doonlys[1].items() == [('task', 'bse')]
