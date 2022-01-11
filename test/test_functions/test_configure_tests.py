"""
Test module that parses YAML config file, and uses it to configure each test setting.

Notes:
 * If one uses a test name containing `properties/test_name`, a tolerance file needs mocking.
   - This is demonstrated in `test_get_method_from_tolerance_file`.
 * Function `configure_all_tests`is not unit-tested because the default input files are found
   from inspecting the test farm => Would also require file mocking.
"""
import pathlib
import pytest
import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from ..src.runner.configure_tests import initialise_tests, assign_test_configurations, get_method_from_test_name, \
    ConfigurationDefaults, get_method_from_tolerance_file, setup_tests, str_attributes_to_enums, \
    parse_config_defaults_file, find_tests_with_attribute, sort_tests_by_group, find_tests_to_skip_by_group

# Renamed else pytest will try and collect from it
from ..src.runner.configure_tests import TestGroup as Group
from ..src.runner.profile import CompilerBuild


def test_str_attributes_enums():
    # GROUP property
    tests = str_attributes_to_enums({'groundstate/LDA_PW-PbTiO3': {'group': 'SLOW_TESTS'}})
    assert tests == {'groundstate/LDA_PW-PbTiO3': {'group': Group.SLOW_TESTS}}

    # Failing builds property
    tests = str_attributes_to_enums({'groundstate/LDA_PW-PbTiO3': {'failing_builds': ['intel_mpiandsmp']}})
    assert tests['groundstate/LDA_PW-PbTiO3']['failing_builds'][0].build == CompilerBuild('intel', 'mpiandsmp').build
    assert tests['groundstate/LDA_PW-PbTiO3']['failing_builds'][0].compiler == CompilerBuild('intel', 'mpiandsmp').compiler

    # Erroneous compiler/build string (missing underscore)
    expected_error = "String representation of CompilerBuild enum is invalid: intel mpiandsmp. \n" \
                     "See src/runner/profile.py or config file documentation for valid strings."

    with pytest.raises(ValueError) as error_info:
        tests = str_attributes_to_enums({'groundstate/LDA_PW-PbTiO3': {'failing_builds': ['intel mpiandsmp']}})

    assert error_info.value.args[0] == expected_error


def test_get_method_from_test_name():
    """
    Test when a method can be extracted from the test case name.
    """

    method = get_method_from_test_name('groundstate/test_name')
    assert method == 'groundstate', "Convention should be method/test_name"

    method = get_method_from_test_name('some/longer/path/groundstate/test_name')
    assert method == 'groundstate', "Longer prepended path"

    method = get_method_from_test_name('properties/test_name')
    assert isinstance(method, KeyError), "properties subdirectory contains multiple methods " \
                               "- cannot infer method from it, hence return KeyError instance"

    method = get_method_from_test_name('junk/test_name')
    assert isinstance(method, KeyError), "Erroneous subdirectory returns KeyError instance"


def test_get_method_from_tolerance_file(tmp_path):
    """
    Test that the method name can be inferred from the tolerance file.
    """

    # Mock directories based on test name
    test_name: pathlib.Path = tmp_path / "si-bandstructure"
    test_name.mkdir()
    ref_directory = test_name / "ref"
    ref_directory.mkdir()
    assert ref_directory.is_dir(), "tmp_path directory does not exist"

    # Mock JSON file with correct name
    file = ref_directory / "tolerance_bandstructure.json"
    file.write_text("{'Arbitrary':'data'")

    assert get_method_from_tolerance_file(test_name.as_posix()) == 'band_structure'


def test_get_method_from_erroneous_tolerance_file(tmp_path):

    test_name: pathlib.Path = tmp_path / "si-bandstructure"
    test_name.mkdir()

    file = test_name / "tolerance_random_method.json"
    file.write_text("{'Arbitrary':'data'")

    # Note, if one applies autoformatting, the assertion on this result will break.
    expected_error = """Invalid tolerance file, tolerance_random_method.json, for test_get_method_from_erroneous0/si-bandstructure. 
Method must correspond to a tolerance:
tolerance_ground_state.json
tolerance_gw.json
tolerance_hybrid.json
tolerance_bse.json
tolerance_xanes.json
tolerance_tddft.json
tolerance_rt_tddft.json
tolerance_dos.json
tolerance_bandstructure.json
tolerance_plotting.json
tolerance_wannier.json
tolerance_transport.json
tolerance_optical.json
tolerance_electric.json
tolerance_core.json
tolerance_spin.json"""

    with pytest.raises(KeyError) as error_info:
        method = get_method_from_tolerance_file(test_name.as_posix(), subdirectory='')

    assert error_info.value.args[0] == expected_error


def test_initialise_tests():
    """
    Test initialisation of the `tests` dictionary, for some choice of default values.
    """
    test_names = ['groundstate/LDA-si', 'gw/LDA-si']
    input_files = ['input.xml', 'si.xml']
    test_defaults = ConfigurationDefaults(files_under_test={'groundstate': ['INFO.OUT'], 'gw': ['GW_INFO.OUT']},
                                          inputs={name: input_files for name in test_names},
                                          repeat=False,
                                          failing_builds=[],
                                          #depends_on=[],
                                          group=Group.NONE,
                                          comments=''
                                          )

    expected_initialisation = \
        {'groundstate/LDA-si':
             {'files_under_test': ['INFO.OUT'],
              'inputs': ['input.xml', 'si.xml'],
              'repeat': False,
              'failing_builds': [],
              #'depends_on': [],
              'group': Group.NONE,
              'comments': ''},
         'gw/LDA-si':
             {'files_under_test': ['GW_INFO.OUT'],
              'inputs': ['input.xml', 'si.xml'],
              'repeat': False,
              'failing_builds': [],
              #'depends_on': [],
              'group': Group.NONE,
              'comments': ''}
         }

    tests = initialise_tests(test_names, test_defaults)
    assert tests == expected_initialisation


def test_initialise_tests_bad_inputs():
    """
    Bad input file key
    """
    test_names = ['groundstate/LDA-si', 'gw/LDA-si']

    input_files_with_bad_key = {'erroneous_method/LDA-si': ['input.xml', 'si.xml'],
                                'gw/LDA-si': ['input.xml', 'si.xml']
                                }
    test_defaults = ConfigurationDefaults(files_under_test={'groundstate': ['INFO.OUT'], 'gw': ['GW_INFO.OUT']},
                                          inputs=input_files_with_bad_key,
                                          repeat=False,
                                          failing_builds=[],
                                          #depends_on=[],
                                          group=Group.NONE,
                                          comments=''
                                          )

    with pytest.raises(KeyError) as error_info:
        tests = initialise_tests(test_names, test_defaults)

    assert error_info.value.args[0] == \
           "Default input keys (from inspecting test_farm) inconsistent with test_names (from config file): " \
           "{'erroneous_method/LDA-si'}"


def test_assign_test_configurations_initialised_are_not_nulled():
    """
     No attributes set in the config file for test.
     Test that the initialised values are not overwritten by null entries.

    Test with all attributes initialised to:
      a) empty/null
      b) Values
    """
    data_from_config = {'groundstate/LDA_PW-PbTiO3': None}

    # Test initialised with null values.
    initialised_tests = {'groundstate/LDA_PW-PbTiO3':
                             {'repeat': False,
                              'failing_builds': [],
                              'files_under_test': [],
                              #'depends_on': '',
                              'group': 'NONE',
                              'inputs': [],
                              'comments': ''
                             }
                         }

    tests = assign_test_configurations(initialised_tests, data_from_config)
    assert tests == initialised_tests, "No attributes specified in `data_from_config` " \
                                       "=> tests should equal the initialised values"

    # Test initialised with values.
    initialised_tests = {'groundstate/LDA_PW-PbTiO3':
                             {'repeat': True,
                              'failing_builds': ['Intel_mpismp'],
                              'files_under_test': ['INFO.OUT'],
                              #'depends_on': None,
                              'group': 'SOME_GROUP',
                              'inputs': ['input.xml', 'Pb.xml', 'Ti.xml', 'O.xml'],
                              'comments': 'Intel MPIANDSMP fails'
                             }
                         }
    tests = assign_test_configurations(initialised_tests, data_from_config)
    assert tests == initialised_tests, "No attributes specified in `data_from_config` " \
                                       "=> tests should equal the initialised values"


@pytest.mark.parametrize("erroneously_initialised_tests", [
                        ({'groundstate/LDA_PW-PbTiO3': {}}),                # All keys missing in initialised test
                        ({'groundstate/LDA_PW-PbTiO3': {'repeat': False}})  # Some keys missing in initialised test
])
def test_assign_test_configurations_erroneous_init(erroneously_initialised_tests):
    """
    Erroneously Initialised Input. Test not initialised with:
     * any attribute keys
     * all attribute keys
    """
    data_from_config = {'groundstate/LDA_PW-PbTiO3': None}

    with pytest.raises(KeyError) as error_info:
        tests = assign_test_configurations(erroneously_initialised_tests, data_from_config)

    assert error_info.value.args[0] == \
           'Test dictionary groundstate/LDA_PW-PbTiO3, not initialised with all ConfigurationDefaults keys'


def test_assign_test_configurations_all_keys_replaced():
    """
    Test expected behaviour. ALl keys in initialised dict replaced by config data.
    """
    initialised_tests = {'groundstate/LDA_PW-PbTiO3':
                             {'repeat': False,
                              'failing_builds': [],
                              'files_under_test': ['INFO.OUT'],
                              #'depends_on': None,
                              'group': 'NULL',
                              'inputs': ['input.xml', 'A.xml', 'B.xml', 'C.xml'],
                              'comments': 'Default comment to be overwritten'
                             }
                         }

    data_from_config = {'groundstate/LDA_PW-PbTiO3':
                             {'repeat': True,
                              'failing_builds': ['Intel_mpismp'],
                              'files_under_test': ['INFO.OUT', 'eigval.xml'],
                              #'depends_on': ['groundstate/another_test'],
                              'group': 'SOME_GROUP',
                              'inputs': ['input.xml', 'Pb.xml', 'Ti.xml', 'O.xml'],
                              'comments': 'Comment specified in config file'
                             }
                         }

    tests = assign_test_configurations(initialised_tests, data_from_config)

    assert tests == data_from_config, "All attributes specified in data_from_config should get added to `tests`" \
                                      "therefore tests should equal the config dictionary"


def test_assign_test_configurations_some_keys_replaced():
    """
    Keys specified in config data replace initialised values.
    Otherwise, initialised values are retained
    """
    initialised_tests = {'groundstate/LDA_PW-PbTiO3':
                             {'repeat': False,
                              'failing_builds': [],
                              'files_under_test': ['INFO.OUT'],
                              #'depends_on': None,
                              'group': 'NULL',
                              'inputs': ['input.xml', 'A.xml', 'B.xml', 'C.xml'],
                              'comments': ''
                             }
                         }

    data_from_config = {'groundstate/LDA_PW-PbTiO3':
                             {'failing_builds': ['Intel_mpismp'],
                              'inputs': ['input.xml', 'Pb.xml', 'Ti.xml', 'O.xml']
                             }
                         }

    tests = assign_test_configurations(initialised_tests, data_from_config)

    expected_tests = {'groundstate/LDA_PW-PbTiO3':
                          {'repeat': False,
                           'failing_builds': ['Intel_mpismp'],
                           'files_under_test': ['INFO.OUT'],
                           #'depends_on': None,
                           'group': 'NULL',
                           'inputs': ['input.xml', 'Pb.xml', 'Ti.xml', 'O.xml'],
                           'comments': ''
                           }
                      }

    assert tests['groundstate/LDA_PW-PbTiO3']['repeat'] == initialised_tests['groundstate/LDA_PW-PbTiO3']['repeat'], \
        "Returned tests should retain initialised values for 'repeat' as it is not specified in config dict"

    assert tests['groundstate/LDA_PW-PbTiO3']['files_under_test'] == \
           initialised_tests['groundstate/LDA_PW-PbTiO3']['files_under_test'], \
        "Returned tests should retain initialised values for 'files_under_test' as it is not specified in config dict"

    # assert tests['groundstate/LDA_PW-PbTiO3']['depends_on'] == \
    #        initialised_tests['groundstate/LDA_PW-PbTiO3']['depends_on'], \
    #     "Returned tests should retain initialised values for 'depends_on' as it is not specified in config dict"

    assert tests['groundstate/LDA_PW-PbTiO3']['group'] == initialised_tests['groundstate/LDA_PW-PbTiO3']['group'], \
        "Returned tests should retain initialised values for 'group' as it is not specified in config dict"

    for key in ['failing_builds', 'inputs']:
        assert tests['groundstate/LDA_PW-PbTiO3'][key] == data_from_config['groundstate/LDA_PW-PbTiO3'][key], \
            "Returned tests should take values when specified in config dict"

    assert tests == expected_tests


@pytest.fixture()
def defaults_config_str() -> str:
    """
    Mock the defaults config file
    """
    defaults_string = """
     default_files_under_test:
        groundstate:
           - "INFO.OUT"
           - "evalcore.xml"
           - "geometry.xml"
           - "eigval.xml"
           - "atoms.xml"
        gw:
           - "some_file"
        hybrid:
           - "some_file"
        bse:
           - "some_file"
        plot:
            - "some_file"
        wannier:
            - "some_file"
        spin_properties:
            - "some_file"
        optical_properties:
            - "some_file"
        tddft:
            - "some_file"
        rt_tddft:
            - "some_file"
        electric_properties:
            - "some_file"
        core_properties:
            - "some_file"
        band_structure:
            - "some_file"
        xanes:
            - "some_file"
        transport:
            - "some_file"
        dos:
            - "some_file"

     # Specify each new group.
     # Any test without a group is assigned NONE
     group_execution:
       NONE: True
       SLOW_TESTS: False
    """
    return defaults_string


@pytest.fixture()
def tests_config_str() -> str:
    """
    Mock the tests config file
    """
    string = """
     # ---------------------------------------------------------------------------------------------------------------
     # exciting Test Suite Configure File
     #
     # This file explicitly specifies each test case.
     #
     # One can explicitly add the attributes/properties:
     #
     #  group:              See defaults file.
     #  repeat:             Repeat if test "sometimes fails". True or False
     #  failing_builds:     For a given test, valid choices are intel_serial, intel_mpiandsmp, gcc_serial, gcc_mpiandsmp
     #  comments:          "Comments regarding any failing tests"
     #  files_under_test:   See defaults file.
     #  inputs:             See defaults file.
     #  depends_on:         Currently not implemented. See issue #
     #
     #  Or if one wishes to use the defaults defined in () simply put the test name:
     #   ```
     #   method_directory/test_name:
     #   ```
     #  where the terminating `:` is required.
     # ---------------------------------------------------------------------------------------------------------------

     # TODO(Alex) Issue 5. Developers can add comments to YAML
     groundstate/LDA_PW-PbTiO3:
        group: SLOW_TESTS
        repeat: False
        failing_builds:
           - intel_mpiandsmp
           - intel_serial
        comments: "Some details on failing test to print from CI"
        files_under_test:
           - "INFO.OUT"
           - "evalcore.xml"
           - "geometry.xml"
           - "eigval.xml"
           - "atoms.xml"

     # Can simply specify the file name, but must include the `:`
     groundstate/LiF:

     # Can specify that a test depends upon another running first.
     # If so, one also must specify which input files to copyk
     BSE/LiF:
        repeat: True
        depends-on: groundstate/LiF
        inputs:
           - input.xml
           - STATE.OUT
           - Li.xml
           - F.xml

        """
    return string


def test_parse_and_setup_tests(defaults_config_str, tests_config_str):
    """
    Take a YAML-formatted string, parse and verify output dictionary
    (without removing any tests)

    Can't easily test `configure_all_tests` because it interacts with the file system
    when finding input files - would need to mock those files.
    """
    config_data: dict = yaml.load(tests_config_str, Loader=Loader)
    config_data = str_attributes_to_enums(config_data)
    config_defaults = parse_config_defaults_file(defaults_config_str)

    mocked_default_inputs = {'groundstate/LDA_PW-PbTiO3': ['input.xml'],
                             'groundstate/LiF': ['input.xml'],
                             'BSE/LiF': ['input.xml'],
                             }

    assert list(config_data.keys()) == list(mocked_default_inputs.keys()), \
        "Test names in `tests_config_str` differ from those specified in `mocked_default_inputs`"

    default_attributes = ConfigurationDefaults(files_under_test=config_defaults['default_files_under_test'],
                                               inputs=mocked_default_inputs,
                                               repeat=False,
                                               failing_builds=[],
                                               #depends_on=[],
                                               group=Group.NONE,
                                               comments=''
                                               )

    tests = setup_tests(config_data, default_attributes)

    expected_tests = {'groundstate/LDA_PW-PbTiO3':
                          {'files_under_test': ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml'],
                           'inputs': ['input.xml'],
                           'repeat': False,
                           'failing_builds': [CompilerBuild('intel', 'mpiandsmp'), CompilerBuild('intel', 'serial')],
                           #'depends_on': [],
                           'group': Group.SLOW_TESTS,
                           'comments': 'Some details on failing test to print from CI'
                           },
                      'groundstate/LiF':
                          {'files_under_test': ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml'],
                           'inputs': ['input.xml'],
                           'repeat': False,
                           'failing_builds': [],
                           #'depends_on': [],
                           'group': Group.NONE,
                           'comments': ''
                           },
                      'BSE/LiF':
                          {'files_under_test': ['some_file'],
                           'inputs': ['input.xml', 'STATE.OUT', 'Li.xml', 'F.xml'],
                           'repeat': True,
                           'failing_builds': [],
                           #'depends_on': [],
                           'group': Group.NONE,
                           'depends-on': 'groundstate/LiF',
                           'comments': ''
                           }
                      }

    # Can't compare different instantiations of an object, even when it has the same attributes set.
    PbTiO3_failing_builds = tests['groundstate/LDA_PW-PbTiO3'].pop('failing_builds')
    expected_PbTiO3_failing_builds = expected_tests['groundstate/LDA_PW-PbTiO3'].pop('failing_builds')

    assert tests == expected_tests

    assert len(PbTiO3_failing_builds) == 2, "Expect two failing build entries for groundstate/LDA_PW-PbTiO3"

    # Entry one
    assert PbTiO3_failing_builds[0].compiler == expected_PbTiO3_failing_builds[0].compiler
    assert PbTiO3_failing_builds[0].build == expected_PbTiO3_failing_builds[0].build

    # Entry two
    assert PbTiO3_failing_builds[1].compiler == expected_PbTiO3_failing_builds[1].compiler
    assert PbTiO3_failing_builds[1].build == expected_PbTiO3_failing_builds[1].build


@pytest.mark.parametrize("key, default, expected_tests", [
                        ('failing_builds', [], ['groundstate/LDA_PW-PbTiO3']),
                        ('repeat', False, ['BSE/LiF'])
])
def test_find_tests_with_attribute(key, default, expected_tests):
    """
    Find tests with attributes that don not equal their default (initialised) values.

    Find test cases the attribute:
        failing_builds
        repeat
    """
    tests_from_config = {
    'groundstate/LDA_PW-PbTiO3': {'repeat': False,
                                  'failing_builds': ['intel mpiandsmp', 'intel serial'],
                                  'files_under_test': ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml','atoms.xml'],
                                  #'depends_on': '',
                                  'group': 'SLOW_TESTS',
                                  'comments': 'Some details on failing test to print from CI'},

     'groundstate/LiF': {'repeat': False,
                         'failing_builds': [],
                         'files_under_test': ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml'],
                         'depends-on': '',
                         'group': 'NONE'},

     'BSE/LiF': {'repeat': True,
                 'failing_builds': [],
                 'files_under_test': ['some_file'],
                 #'depends_on': 'groundstate/LiF',
                 'group': 'NONE',
                 'inputs': ['input.xml', 'STATE.OUT', 'Li.xml', 'F.xml']
                 }
    }
    output_tests: dict = find_tests_with_attribute(tests_from_config, key, default)
    assert list(output_tests.keys()) == expected_tests


def test_sort_tests_by_group():

    yaml_specified_tests = {
     'groundstate/LDA_PW-PbTiO3': {'group': Group.SLOW_TESTS},
     'groundstate/LiF': {'group': Group.NONE},
     'BSE/LiF': {'group': Group.NONE}
    }

    tests_by_group = sort_tests_by_group(yaml_specified_tests)

    expected_test_grouping = {Group.NONE: ['groundstate/LiF', 'BSE/LiF'],
                              Group.SLOW_TESTS: ['groundstate/LDA_PW-PbTiO3']
                              }

    assert tests_by_group == expected_test_grouping


def test_find_tests_to_skip_by_group():
    """
    List tests with group execution set to false
    """

    yaml_specified_tests = {
     'groundstate/LDA_PW-PbTiO3': {'group': Group.SLOW_TESTS},
     'groundstate/LiF': {'group': Group.NONE},
     'BSE/LiF': {'group': Group.NONE}
    }

    # List SLOW_TESTS
    group_execution = {'NONE': True, 'SLOW_TESTS': False}
    tests_to_skip = find_tests_to_skip_by_group(yaml_specified_tests, group_execution)
    assert tests_to_skip == {"groundstate/LDA_PW-PbTiO3"}, "Expect to list tests tagged with group != None"

    # List NONE
    group_execution = {'NONE': False, 'SLOW_TESTS': True}
    tests_to_skip = find_tests_to_skip_by_group(yaml_specified_tests, group_execution)
    assert tests_to_skip == {'groundstate/LiF', 'BSE/LiF'}, "Expect to list tests tagged with group != SLOW_TESTS"

    # All tests skipped
    group_execution = {'NONE': False, 'SLOW_TESTS': False}
    tests_to_skip = find_tests_to_skip_by_group(yaml_specified_tests, group_execution)
    assert tests_to_skip == {"groundstate/LDA_PW-PbTiO3", 'groundstate/LiF', 'BSE/LiF'}, "Expect all tests to be listed"
