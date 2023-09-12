"""Test ExcitingGroundStateInput class attribute assignment and methods.

NOTE:
All attribute tests should assert on the XML tree content's as the attribute
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
import pytest

from excitingtools.input.input_classes import ExcitingGroundStateInput


@pytest.mark.parametrize(
    "test_input,expected",
    [({"rgkmax": 8.6}, [('rgkmax', '8.6')]), ({"ngridk": [8, 8, 8]}, [('ngridk', '8 8 8')]),
     ({'vkloff': [0.1, 0.2, 0.3]}, [('vkloff', '0.1 0.2 0.3')]),
     ({'CoreRelativity': 'dirac'}, [('CoreRelativity', 'dirac')]),
     ({'ExplicitKineticEnergy': True}, [('ExplicitKineticEnergy', 'true')]),
     ({'PrelimLinSteps': 2}, [('PrelimLinSteps', '2')]),
     ({'ValenceRelativity': 'zora'}, [('ValenceRelativity', 'zora')]),
     ({'autokpt': False}, [('autokpt', 'false')]), ({'beta0': 0.4}, [('beta0', '0.4')]),
     ({'betadec': 0.6}, [('betadec', '0.6')]), ({'betainc': 1.1}, [('betainc', '1.1')]),
     ({'cfdamp': 0.0}, [('cfdamp', '0.0')]), ({'chgexs': 0.0}, [('chgexs', '0.0')]),
     ({'deband': 2.5e-3}, [('deband', '0.0025')]),
     ({'dipolecorrection': True}, [('dipolecorrection', 'true')]),
     ({'dipoleposition': 1.0}, [('dipoleposition', '1.0')]),
     ({'dlinengyfermi': -0.1}, [('dlinengyfermi', '-0.1')]),
     ({'do': "fromscratch"}, [('do', "fromscratch")]),
     ({'energyref': 0.0}, [('energyref', '0.0')]), ({'epsband': 1.0e-6}, [('epsband', '1e-06')]),
     ({'epschg': 1.0e-5}, [('epschg', '1e-05')]), ({'epsengy': 1e-6}, [('epsengy', '1e-06')]),
     ({'epsforcescf': 5.0e-5}, [('epsforcescf', '5e-05')]),
     ({'epsocc': 1e-8}, [('epsocc', '1e-08')]), ({'epspot': 1e-6}, [('epspot', '1e-06')]),
     ({'fermilinengy': False}, [('fermilinengy', 'false')]),
     ({'findlinentype': "Wigner_Seitz"}, [('findlinentype', "Wigner_Seitz")]),
     ({'fracinr': 0.02}, [('fracinr', '0.02')]), ({'frozencore': False}, [('frozencore', 'false')]),
     ({'gmaxvr': 12}, [('gmaxvr', '12')]), ({'isgkmax': -1}, [('isgkmax', '-1')]),
     ({'ldapu': "none"}, [('ldapu', "none")]), ({'lmaxapw': 8}, [('lmaxapw', '8')]),
     ({'lmaxinr': 2}, [('lmaxinr', '2')]), ({'lmaxmat': 8}, [('lmaxmat', '8')]),
     ({'lmaxvr': 8}, [('lmaxvr', '8')]), ({'lradstep': 1}, [('lradstep', '1')]),
     ({'maxscl': 200}, [('maxscl', '200')]), ({'mixer': 'msec'}, [('mixer', 'msec')]),
     ({'mixerswitch': 1}, [('mixerswitch', '1')]), ({'modifiedsv': False}, [('modifiedsv', 'false')]),
     ({'msecStoredSteps': 8}, [('msecStoredSteps', '8')]), ({'nempty': 5}, [('nempty', '5')]),
     ({'niterconvcheck': 2}, [('niterconvcheck', '2')]), ({'nktot': 0}, [('nktot', '0')]),
     ({'nosource': False}, [('nosource', 'false')]), ({'nosym': False}, [('nosym', 'false')]),
     ({'nprad': 4}, [('nprad', '4')]), ({'npsden': 9}, [('npsden', '9')]),
     ({'nwrite': 0}, [('nwrite', '0')]), ({'outputlevel': 'normal'}, [('outputlevel', 'normal')]),
     ({'ptnucl': True}, [('ptnucl', 'true')]), ({'radialgridtype': 'cubic'}, [('radialgridtype', 'cubic')]),
     ({'radkpt': 40.0}, [('radkpt', '40.0')]), ({'reducek': True}, [('reducek', 'true')]),
     ({'scfconv': 'multiple'}, [('scfconv', 'multiple')]), ({'stype': 'Gaussian'}, [('stype', 'Gaussian')]),
     ({'swidth': 0.001}, [('swidth', '0.001')]), ({'symmorph': False}, [('symmorph', 'false')]),
     ({'tevecsv': False}, [('tevecsv', 'false')]), ({'tfibs': True}, [('tfibs', 'true')]),
     ({'tforce': False}, [('tforce', 'false')]), ({'tpartcharges': False}, [('tpartcharges', 'false')]),
     ({'useDensityMatrix': True}, [('useDensityMatrix', 'true')]),
     ({'vdWcorrection': 'none'}, [('vdWcorrection', 'none')]),
     ({'xctype': 'GGA_PBE'}, [('xctype', 'GGA_PBE')])]
)
def test_class_exciting_ground_state_input_parametrized(test_input, expected):
    gs_input = ExcitingGroundStateInput(**test_input)
    gs_xml = gs_input.to_xml()
    assert gs_xml.tag == 'groundstate'
    assert gs_xml.items() == expected


def test_invalid_input():
    """
    Test error is raised when giving bogus attributes to class constructor.
    """
    # Use an erroneous ground state attribute
    with pytest.raises(ValueError, match="groundstate keys are not valid: {'erroneous_attribute'}"):
        ExcitingGroundStateInput(erroneous_attribute=True)


def test_as_dict(mock_env_jobflow_missing):
    ref_rgkmax = 8.5
    gs_input = ExcitingGroundStateInput(rgkmax=ref_rgkmax)
    ref_dict = {'xml_string': f'<groundstate rgkmax="{ref_rgkmax}"> </groundstate>'}
    assert gs_input.as_dict() == ref_dict, 'expected different dict representation'


def test_as_dict_jobflow(mock_env_jobflow):
    ref_rgkmax = 8.5
    gs_input = ExcitingGroundStateInput(rgkmax=ref_rgkmax)
    ref_dict = {'@class': 'ExcitingGroundStateInput',
                '@module': 'excitingtools.input.input_classes',
                'xml_string': f'<groundstate rgkmax="{ref_rgkmax}"> </groundstate>'}
    assert gs_input.as_dict() == ref_dict, 'expected different dict representation'


def test_from_dict():
    ref_rgkmax = 8.5
    ref_dict = {'xml_string': f'<groundstate rgkmax="{ref_rgkmax}"> </groundstate>'}
    gs_input = ExcitingGroundStateInput.from_dict(ref_dict)

    assert gs_input.name == 'groundstate'
    # added comment for pylint to disable warning, because of dynamic attributes
    assert gs_input.rgkmax == ref_rgkmax, f'Expect rgkmax to be equal {ref_rgkmax}'  # pylint: disable=no-member


def test_spin_input():
    spin_attributes = {'bfieldc': [0, 0, 0], 'fixspin': 'total FSM'}
    spin_keys = list(spin_attributes)
    gs_input = ExcitingGroundStateInput(rgkmax=7.0, spin=spin_attributes)

    gs_xml = gs_input.to_xml()
    assert gs_xml.tag == 'groundstate'
    assert set(gs_xml.keys()) == {'rgkmax'}

    elements = list(gs_xml)
    assert len(elements) == 1

    spin_xml = elements[0]
    assert spin_xml.tag == "spin"
    assert spin_xml.keys() == spin_keys, 'Should contain all spin attributes'
    assert spin_xml.get('bfieldc') == '0 0 0'
    assert spin_xml.get('fixspin') == 'total FSM'


def test_solver_input():
    solver_attributes = {'packedmatrixstorage': True, 'type': "Lapack"}
    solver_keys = list(solver_attributes)
    gs_input = ExcitingGroundStateInput(solver=solver_attributes)

    gs_xml = gs_input.to_xml()
    assert gs_xml.tag == 'groundstate'
    assert set(gs_xml.keys()) == set()

    elements = list(gs_xml)
    assert len(elements) == 1

    solver_xml = elements[0]
    assert solver_xml.tag == "solver"
    assert solver_xml.keys() == solver_keys, 'Should contain all spin attributes'
    assert solver_xml.get('packedmatrixstorage') == 'true'
    assert solver_xml.get('type') == 'Lapack'


def test_dfthalf_input():
    dfthalf_attributes = {"printVSfile": True}
    dfthalf_keys = list(dfthalf_attributes)
    gs_input = ExcitingGroundStateInput(dfthalf=dfthalf_attributes)

    gs_xml = gs_input.to_xml()
    assert gs_xml.tag == 'groundstate'
    assert set(gs_xml.keys()) == set()

    elements = list(gs_xml)
    assert len(elements) == 1

    dfthalf_xml = elements[0]
    assert dfthalf_xml.tag == "dfthalf"
    assert dfthalf_xml.keys() == dfthalf_keys, 'Should contain all dfthalf attributes'
    assert dfthalf_xml.get('printVSfile') == 'true'


def test_OEP_input():
    oep_attributes = {"convoep": 1e-5, "maxitoep": 200, "tauoep": [1, 2, 3]}
    oep_keys = list(oep_attributes)
    gs_input = ExcitingGroundStateInput(OEP=oep_attributes)

    gs_xml = gs_input.to_xml()
    assert gs_xml.tag == 'groundstate'
    assert set(gs_xml.keys()) == set()

    elements = list(gs_xml)
    assert len(elements) == 1

    oep_xml = elements[0]
    assert oep_xml.tag == "OEP"
    assert oep_xml.keys() == oep_keys, 'Should contain all dfthalf attributes'
    assert oep_xml.get('convoep') == '1e-05'
    assert oep_xml.get('maxitoep') == '200'
    assert oep_xml.get('tauoep') == '1 2 3'
