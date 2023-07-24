"""Test ExcitingGroundStateInput class attribute assignment and methods.

NOTE:
All attribute tests should assert on the XML tree content's as the attribute
order is not preserved by the ElementTree.tostring method. Elements appear to
be fine.
"""
import numpy as np
import pytest

# noinspection PyUnresolvedReferences
from excitingtools.input.input_classes import ExcitingBSEInput, ExcitingXSInput


def test_class_ExcitingXSInput():
    xs_input = ExcitingXSInput(xstype="BSE")

    xs_xml = xs_input.to_xml()
    assert xs_xml.tag == 'xs'
    assert xs_xml.items() == [('xstype', 'BSE')]


def test_class_ExcitingXSInput_xstype_missing():
    with pytest.raises(ValueError, match="Missing mandatory arguments: {'xstype'}"):
        ExcitingXSInput(ngridk=[4, 4, 4])


def test_ExcitingXSInput_as_dict(mock_env_jobflow_missing):
    xs_input = ExcitingXSInput(xstype="BSE")
    ref_dict = {'xml_string': '<xs xstype="BSE"> </xs>'}
    assert xs_input.as_dict() == ref_dict, 'expected different dict representation'


def test_ExcitingXSInput_as_dict_jobflow(mock_env_jobflow):
    xs_input = ExcitingXSInput(xstype="BSE")
    ref_dict = {'@class': 'ExcitingXSInput',
                '@module': 'excitingtools.input.input_classes',
                'xml_string': '<xs xstype="BSE"> </xs>'}
    assert xs_input.as_dict() == ref_dict, 'expected different dict representation'


def test_ExcitingXSInput_from_dict():
    ref_dict = {'xml_string': '<xs xstype="BSE"> </xs>'}
    recreated_xs = ExcitingXSInput.from_dict(ref_dict)
    assert recreated_xs.to_xml_str() == ref_dict['xml_string']


def test_class_ExcitingXSInput_xs():
    xs = {'broad': 0.32, 'ngridk': [8, 8, 8], 'tevout': True, 'nempty': 52, 'pwmat': 'fft', "xstype": "BSE"}
    xs_input = ExcitingXSInput(**xs)

    xs_xml = xs_input.to_xml()
    assert xs_xml.tag == 'xs'
    assert set(xs_xml.keys()) == set(xs)
    assert xs_xml.get('xstype') == 'BSE'
    assert xs_xml.get('broad') == '0.32'
    assert xs_xml.get('ngridk') == '8 8 8'
    assert xs_xml.get('tevout') == 'true'
    assert xs_xml.get('nempty') == '52'
    assert xs_xml.get('pwmat') == 'fft'


def test_class_ExcitingXsInput_wrong_key():
    with pytest.raises(ValueError, match="xs keys are not valid: {'wrong_key'}"):
        ExcitingXSInput(xstype="BSE", wrong_key=1)


def test_class_ExcitingXSInput_BSE_element():
    bse_attributes = {'bsetype': 'singlet', 'xas': True, 'xasspecies': 1}
    bse_keys = list(bse_attributes)
    xs_input = ExcitingXSInput(xstype="BSE", BSE=bse_attributes)

    xs_xml = xs_input.to_xml()
    assert xs_xml.tag == 'xs'
    assert set(xs_xml.keys()) == {'xstype'}

    elements = list(xs_xml)
    assert len(elements) == 1

    bse_xml = elements[0]
    assert bse_xml.tag == "BSE"
    assert bse_xml.keys() == bse_keys, 'Should contain all bse attributes'
    assert bse_xml.get('bsetype') == 'singlet'
    assert bse_xml.get('xas') == 'true'
    assert bse_xml.get('xasspecies') == '1'


def test_class_ExcitingXSInput_BSE_element_object():
    bse_object = ExcitingBSEInput(bsetype='singlet', xas=True, xasspecies=1)
    bse_keys = {'bsetype', 'xas', 'xasspecies'}
    xs_input = ExcitingXSInput(xstype="BSE", BSE=bse_object)

    xs_xml = xs_input.to_xml()
    assert xs_xml.tag == 'xs'
    assert set(xs_xml.keys()) == {'xstype'}

    elements = list(xs_xml)
    assert len(elements) == 1

    bse_xml = elements[0]
    assert bse_xml.tag == "BSE"
    assert set(bse_xml.keys()) == bse_keys, 'Should contain all bse attributes'
    assert bse_xml.get('bsetype') == 'singlet'
    assert bse_xml.get('xas') == 'true'
    assert bse_xml.get('xasspecies') == '1'


def test_class_ExcitingXSInput_energywindow_element():
    energywindow_attributes = {'intv': [5.8, 8.3], 'points': 5000}
    xs_input = ExcitingXSInput(xstype="BSE", energywindow=energywindow_attributes)

    xs_xml = xs_input.to_xml()
    assert xs_xml.tag == 'xs'

    elements = list(xs_xml)
    assert len(elements) == 1

    energywindow_xml = elements[0]
    assert energywindow_xml.tag == "energywindow"
    assert set(energywindow_xml.keys()) == set(energywindow_attributes), 'Should contain all bse attributes'
    assert energywindow_xml.get('intv') == '5.8 8.3'
    assert energywindow_xml.get('points') == '5000'


def test_class_ExcitingXSInput_screening_element():
    screening_attributes = {'screentype': 'full', 'nempty': 15}
    xs_input = ExcitingXSInput(xstype="BSE", screening=screening_attributes)

    xs_xml = xs_input.to_xml()
    assert xs_xml.tag == 'xs'

    elements = list(xs_xml)
    assert len(elements) == 1

    screening_xml = elements[0]
    assert screening_xml.tag == "screening"
    assert set(screening_xml.keys()) == set(screening_attributes), 'Should contain all bse attributes'
    assert screening_xml.get('screentype') == 'full'
    assert screening_xml.get('nempty') == '15'


def test_class_ExcitingQpointsetInput_numpy():
    qpointset_input = np.array(((0, 0, 0), (0.5, 0.5, 0.5)))

    xs_input = ExcitingXSInput(xstype="BSE", qpointset=qpointset_input)

    xs_xml = xs_input.to_xml()
    assert xs_xml.tag == 'xs'

    elements = list(xs_xml)
    assert len(elements) == 1

    qpointset_xml = elements[0]
    assert qpointset_xml.tag == "qpointset"
    assert qpointset_xml.items() == []

    qpoints = list(qpointset_xml)
    assert len(qpoints) == 2
    assert qpoints[0].tag == 'qpoint'
    assert qpoints[0].items() == []
    assert qpoints[0].text == '0.0 0.0 0.0'
    assert qpoints[1].text == '0.5 0.5 0.5'


def test_class_ExcitingQpointsetInput_list():
    qpointset_input = [[0, 0, 0], [0.5, 0.5, 0.5]]

    xs_input = ExcitingXSInput(xstype="BSE", qpointset=qpointset_input)

    xs_xml = xs_input.to_xml()
    assert xs_xml.tag == 'xs'

    elements = list(xs_xml)
    assert len(elements) == 1

    qpointset_xml = elements[0]
    assert qpointset_xml.tag == "qpointset"
    assert qpointset_xml.items() == []

    qpoints = list(qpointset_xml)
    assert len(qpoints) == 2
    assert qpoints[0].tag == 'qpoint'
    assert qpoints[0].items() == []
    assert qpoints[0].text == '0 0 0'
    assert qpoints[1].text == '0.5 0.5 0.5'


def test_class_ExcitingPlanInput():
    plan_input = ['screen', 'bse', 'bsegenspec']

    xs_input = ExcitingXSInput(xstype="BSE", plan=plan_input)

    xs_xml = xs_input.to_xml()
    assert xs_xml.tag == 'xs'

    elements = list(xs_xml)
    assert len(elements) == 1

    plan_xml = elements[0]
    assert plan_xml.tag == "plan"
    assert plan_xml.items() == []

    doonlys = list(plan_xml)
    assert len(doonlys) == 3
    assert doonlys[0].tag == 'doonly'
    assert doonlys[0].items() == [('task', 'screen')]
    assert doonlys[1].tag == 'doonly'
    assert doonlys[1].items() == [('task', 'bse')]
    assert doonlys[2].tag == 'doonly'
    assert doonlys[2].items() == [('task', 'bsegenspec')]


def test_class_ExcitingPlanInput_wrong_plan():
    plan_input = ['screen', 'bse', 'bsegenspec', 'falseplan']

    with pytest.raises(ValueError) as error:
        ExcitingXSInput(xstype="BSE", plan=plan_input)
    assert error.value.args[0] == "plan keys are not valid: {'falseplan'}"


def test_class_ExcitingXSInput_attribute_setting_getting():
    xs_input = ExcitingXSInput(xstype="BSE")

    xs_input.ngridk = [2, 2, 2]
    with pytest.raises(ValueError, match="xs keys are not valid: {'abc'}"):
        xs_input.abc = 3


def test_class_ExcitingXSInput_attribute_deleting():
    xs_input = ExcitingXSInput(xstype="BSE", ngridk=[2, 2, 2])

    del xs_input.ngridk
    with pytest.warns(UserWarning, match="Attempt to delete mandatory attribute 'xstype' was prevented."):
        del xs_input.xstype  # pylint: disable=no-member
