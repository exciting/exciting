"""Test ExcitingGroundStateInput class attribute assignment and methods
"""
import numpy as np
import pytest

from excitingtools.input.xs import ExcitingXSInput
from excitingtools.input.base_class import ExcitingXMLInput


def test_class_ExcitingXSInput():
    xs_input = ExcitingXSInput("BSE")

    xs_xml = xs_input.to_xml()
    assert xs_xml.tag == 'xs'
    assert xs_xml.items() == [('xstype', 'BSE')]


def test_class_ExcitingXSInput_xs():
    xs = {'broad': 0.32, 'ngridk': [8, 8, 8], 'tevout': True, 'nempty': 52, 'pwmat': 'fft'}
    mandatory = ['xstype']
    optional = list(xs)
    xs_input = ExcitingXSInput("BSE", xs=xs)

    xs_xml = xs_input.to_xml()
    assert xs_xml.tag == 'xs'
    try:
        assert xs_xml.keys() == mandatory + optional
    except AssertionError:
        assert xs_xml.keys() == optional + mandatory, 'Should contain mandatory xstype plus all optional attributes'
    assert xs_xml.get('xstype') == 'BSE'
    assert xs_xml.get('broad') == '0.32'
    assert xs_xml.get('ngridk') == '8 8 8'
    assert xs_xml.get('tevout') == 'true'
    assert xs_xml.get('nempty') == '52'
    assert xs_xml.get('pwmat') == 'fft'


def test_class_ExcitingXsInput_wrong_key():
    with pytest.raises(ValueError) as error:
        ExcitingXSInput("BSE", xs={'wrong_key': 1})
    assert error.value.args[0] == "xs keys are not valid: {'wrong_key'}"


def test_class_ExcitingXSInput_BSE_element():
    bse_attributes = {'bsetype': 'singlet', 'xas': True, 'xasspecies': 1}
    bse_keys = list(bse_attributes)
    xs_input = ExcitingXSInput("BSE", BSE=bse_attributes)

    xs_xml = xs_input.to_xml()
    assert xs_xml.tag == 'xs'
    assert xs_xml.keys() == ['xstype']

    elements = list(xs_xml)
    assert len(elements) == 1

    bse_xml = elements[0]
    assert bse_xml.tag == "BSE"
    assert bse_xml.keys() == bse_keys, 'Should contain all bse attributes'
    assert bse_xml.get('bsetype') == 'singlet'
    assert bse_xml.get('xas') == 'true'
    assert bse_xml.get('xasspecies') == '1'


def test_class_ExcitingXSInput_BSE_element_object():
    bse_object = ExcitingXMLInput('BSE', bsetype='singlet', xas=True, xasspecies=1)
    bse_keys = ['bsetype', 'xas', 'xasspecies']
    xs_input = ExcitingXSInput("BSE", BSE=bse_object)

    xs_xml = xs_input.to_xml()
    assert xs_xml.tag == 'xs'
    assert xs_xml.keys() == ['xstype']

    elements = list(xs_xml)
    assert len(elements) == 1

    bse_xml = elements[0]
    assert bse_xml.tag == "BSE"
    assert bse_xml.keys() == bse_keys, 'Should contain all bse attributes'
    assert bse_xml.get('bsetype') == 'singlet'
    assert bse_xml.get('xas') == 'true'
    assert bse_xml.get('xasspecies') == '1'


def test_class_ExcitingXSInput_energywindow_element():
    energywindow_attributes = {'intv': [5.8, 8.3], 'points': 5000}
    energywindow_keys = list(energywindow_attributes)
    xs_input = ExcitingXSInput("BSE", energywindow=energywindow_attributes)

    xs_xml = xs_input.to_xml()
    assert xs_xml.tag == 'xs'

    elements = list(xs_xml)
    assert len(elements) == 1

    energywindow_xml = elements[0]
    assert energywindow_xml.tag == "energywindow"
    assert energywindow_xml.keys() == energywindow_keys, 'Should contain all bse attributes'
    assert energywindow_xml.get('intv') == '5.8 8.3'
    assert energywindow_xml.get('points') == '5000'


def test_class_ExcitingXSInput_screening_element():
    screening_attributes = {'screentype': 'full', 'nempty': 15}
    screening_keys = list(screening_attributes)
    xs_input = ExcitingXSInput("BSE", screening=screening_attributes)

    xs_xml = xs_input.to_xml()
    assert xs_xml.tag == 'xs'

    elements = list(xs_xml)
    assert len(elements) == 1

    screening_xml = elements[0]
    assert screening_xml.tag == "screening"
    assert screening_xml.keys() == screening_keys, 'Should contain all bse attributes'
    assert screening_xml.get('screentype') == 'full'
    assert screening_xml.get('nempty') == '15'


def test_class_ExcitingQpointsetInput_numpy():
    qpointset_input = np.array(((0, 0, 0), (0.5, 0.5, 0.5)))

    xs_input = ExcitingXSInput("BSE", qpointset=qpointset_input)

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

    xs_input = ExcitingXSInput("BSE", qpointset=qpointset_input)

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

    xs_input = ExcitingXSInput("BSE", plan=plan_input)

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
        ExcitingXSInput("BSE", plan=plan_input)
    assert error.value.args[0] == "Plan keys are not valid: {'falseplan'}"
