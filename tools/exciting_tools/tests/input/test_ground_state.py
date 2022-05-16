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

from excitingtools.input.ground_state import ExcitingGroundStateInput


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
     ({'deband': 2.5e-3}, [('deband', '0.0025')]), ({'dipolecorrection': True}, [
         ('dipolecorrection', 'true')
         ]), ({'dipoleposition': 1.0}, [('dipoleposition', '1.0')]), ({'dlinengyfermi': -0.1}, [
             ('dlinengyfermi', '-0.1')
             ]), ({'do': "fromscratch"}, [('do', "fromscratch")]),
     ({'energyref': 0.0}, [('energyref', '0.0')]), ({'epsband': 1.0e-6}, [('epsband', '1e-06')]),
     ({'epschg': 1.0e-5}, [('epschg', '1e-05')]), ({'epsengy': 1e-6}, [('epsengy', '1e-06')]),
     ({'epsforcescf': 5.0e-5}, [('epsforcescf', '5e-05')]),
     ({'epsocc': 1e-8}, [('epsocc', '1e-08')]), ({'epspot': 1e-6}, [('epspot', '1e-06')]),
     ({'fermilinengy': False}, [('fermilinengy', 'false')]), ({'findlinentype': "Wigner_Seitz"}, [
         ('findlinentype', "Wigner_Seitz")
         ]), ({'fracinr': 0.02}, [('fracinr', '0.02')]), ({'frozencore': False}, [
             ('frozencore', 'false')
             ]), ({'gmaxvr': 12}, [('gmaxvr', '12')]), ({'isgkmax': -1}, [('isgkmax', '-1')]),
     ({'ldapu': "none"}, [('ldapu', "none")]), ({'lmaxapw': 8}, [('lmaxapw', '8')]),
     ({'lmaxinr': 2}, [('lmaxinr', '2')]), ({'lmaxmat': 8}, [('lmaxmat', '8')]), ({'lmaxvr': 8}, [
         ('lmaxvr', '8')
         ]), ({'lorecommendation': False}, [('lorecommendation', 'false')]), ({'lradstep': 1}, [
             ('lradstep', '1')
             ]), ({'maxscl': 200}, [('maxscl', '200')]), ({'mixer': 'msec'}, [('mixer', 'msec')]),
     ({'mixerswitch': 1}, [('mixerswitch', '1')]), ({'modifiedsv': False}, [
         ('modifiedsv', 'false')
         ]), ({'msecStoredSteps': 8}, [('msecStoredSteps', '8')]), ({'nempty': 5}, [
             ('nempty', '5')
             ]), ({'niterconvcheck': 2}, [('niterconvcheck', '2')]), ({'nktot': 0}, [
                 ('nktot', '0')
                 ]), ({'nosource': False}, [('nosource', 'false')]), ({'nosym': False
                                                                       }, [('nosym', 'false')]),
     ({'nprad': 4}, [('nprad', '4')]), ({'npsden': 9}, [('npsden', '9')]), ({'nwrite': 0}, [
         ('nwrite', '0')
         ]), ({'outputlevel': 'normal'}, [('outputlevel', 'normal')]), ({'ptnucl': True}, [
             ('ptnucl', 'true')
             ]), ({'radialgridtype': 'cubic'}, [('radialgridtype', 'cubic')]), ({'radkpt': 40.0}, [
                 ('radkpt', '40.0')
                 ]), ({'reducek': True}, [('reducek', 'true')]), ({'scfconv': 'multiple'}, [
                     ('scfconv', 'multiple')
                     ]), ({'stype': 'Gaussian'}, [('stype', 'Gaussian')]), ({'swidth': 0.001}, [
                         ('swidth', '0.001')
                         ]), ({'symmorph': False}, [('symmorph', 'false')]), ({'tevecsv': False}, [
                             ('tevecsv', 'false')
                             ]), ({'tfibs': True}, [('tfibs', 'true')]), ({'tforce': False}, [
                                 ('tforce', 'false')
                                 ]), ({'tpartcharges': False}, [('tpartcharges', 'false')]),
     ({'useDensityMatrix': True}, [('useDensityMatrix', 'true')]), ({'vdWcorrection': 'none'}, [
         ('vdWcorrection', 'none')
         ]), ({'xctype': 'GGA_PBE'}, [('xctype', 'GGA_PBE')])]
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
    with pytest.raises(ValueError) as error:
        ExcitingGroundStateInput(erroneous_attribute=True)

    assert error.value.args[
               0] == "groundstate keys are not valid: {'erroneous_attribute'}"
