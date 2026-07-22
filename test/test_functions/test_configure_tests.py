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

from ..src.io.yaml_configuration import Group, yaml_load, group_tag, enum_group_constructor, build_tag, \
    enum_compiler_build_constructor
from ..src.runner.profile import CompilerBuild, Compiler, BuildType

from ..src.runner.configure_tests import initialise_tests, get_method_from_test_name, get_method_from_tolerance_file, \
    assign_config_keys, find_tests_with_attribute, sort_tests_by_group, find_tests_to_skip_by_group


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
tolerance_bse_hdf5.json
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
tolerance_spin.json
tolerance_wfplot.json
tolerance_fastBSE.json
tolerance_phonon.json"""

    with pytest.raises(KeyError) as error_info:
        method = get_method_from_tolerance_file(test_name.as_posix(), subdirectory='')

    assert error_info.value.args[0] == expected_error


def test_initialise_tests():
    """
    Test initialisation of the `tests` dictionary, for some choice of default values.
    Note, most defaults are defined in `io.yaml_configuration` namedTupled.
    """
    test_names = ['groundstate/LDA-si', 'gw/LDA-si']
    input_files = ['input.xml', 'si.xml']

    default_inputs = {name: input_files for name in test_names}
    default_files_under_test = {'groundstate': ['INFO.OUT'], 'gw': ['GW_INFO.OUT']}

    expected_initialisation = \
        {'groundstate/LDA-si':
             {'files_under_test': ['INFO.OUT'],
              'inputs': ['input.xml', 'si.xml'],
              'repeat': False,
              'failing_builds': [],
              'depends_on': None,
              'group': Group.NONE,
              'comments': '',
              'cmd_line_args': ''},
         'gw/LDA-si':
             {'files_under_test': ['GW_INFO.OUT'],
              'inputs': ['input.xml', 'si.xml'],
              'repeat': False,
              'failing_builds': [],
              'depends_on': None,
              'group': Group.NONE,
              'comments': '',
              'cmd_line_args': ''}
         }

    tests = initialise_tests(test_names, default_inputs, default_files_under_test)
    assert tests == expected_initialisation


def test_initialise_tests_bad_inputs():
    """
    Bad input file key
    """
    test_names = ['groundstate/LDA-si', 'gw/LDA-si']

    input_files_with_bad_key = {'erroneous_method/LDA-si': ['input.xml', 'si.xml'],
                                'gw/LDA-si': ['input.xml', 'si.xml']
                                }
    default_files_under_test = {'groundstate': ['INFO.OUT'], 'gw': ['GW_INFO.OUT']}

    with pytest.raises(KeyError) as error_info:
        tests = initialise_tests(test_names, input_files_with_bad_key, default_files_under_test)

    assert error_info.value.args[0] == \
           "Default input keys (from inspecting test_farm) inconsistent with test_names (from config file): " \
           "{'erroneous_method/LDA-si'}"


@pytest.fixture()
def defaults_config_dict() -> dict:
    """
    Mock the defaults config file
    """
    # Do not include all methods, and add dummy output file names for some
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

     # Specify each new group.
     # Any test without a group is assigned NONE
     group_execution:
       NONE: True
       SLOW_TESTS: False
    """
    return yaml_load(defaults_string)


@pytest.fixture()
def tests_config_dict() -> dict:
    """
    Mock the tests config file, and test its parsing
    """
    string = """
     # exciting Test Suite Configure File

     # TODO(Alex) Issue 5. Developers can add comments to YAML
     groundstate/LDA_PW-PbTiO3:
        group: !Group SLOW_TESTS
        repeat: False
        failing_builds:
           - !Build intel_mpiandsmp
           - !Build intel_serial
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
        depends_on: groundstate/LiF
        inputs:
           - input.xml
           - STATE.OUT
           - Li.xml
           - F.xml

        """
    custom_constructor = {group_tag: enum_group_constructor, build_tag: enum_compiler_build_constructor}
    config = yaml_load(string, custom_constructor=custom_constructor)

    # Test parsing, with enum conversion
    expected_config = {'groundstate/LDA_PW-PbTiO3':
                           {'group': Group.SLOW_TESTS,
                            'repeat': False,
                            'failing_builds': [CompilerBuild('intel', 'mpiandsmp'), CompilerBuild('intel', 'serial')],
                            'comments': 'Some details on failing test to print from CI',
                            'files_under_test': ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml']
                            },
                       'groundstate/LiF': None,
                       'BSE/LiF':
                           {'repeat': True,
                            'depends_on': 'groundstate/LiF',
                            'inputs': ['input.xml', 'STATE.OUT', 'Li.xml', 'F.xml']
                            }
                       }

    # Test name 1
    assert config['groundstate/LDA_PW-PbTiO3']['group'] == expected_config['groundstate/LDA_PW-PbTiO3']['group']
    assert config['groundstate/LDA_PW-PbTiO3']['repeat'] == expected_config['groundstate/LDA_PW-PbTiO3']['repeat']
    assert config['groundstate/LDA_PW-PbTiO3']['comments'] == expected_config['groundstate/LDA_PW-PbTiO3']['comments']
    assert config['groundstate/LDA_PW-PbTiO3']['files_under_test'] == expected_config['groundstate/LDA_PW-PbTiO3']['files_under_test']
    # Cannot directly compare object instantiations
    assert len(config['groundstate/LDA_PW-PbTiO3']['failing_builds']) == 2
    assert config['groundstate/LDA_PW-PbTiO3']['failing_builds'][0].compiler == Compiler.intel
    assert config['groundstate/LDA_PW-PbTiO3']['failing_builds'][0].build == BuildType.mpiandsmp
    assert config['groundstate/LDA_PW-PbTiO3']['failing_builds'][1].compiler == Compiler.intel
    assert config['groundstate/LDA_PW-PbTiO3']['failing_builds'][1].build == BuildType.serial

    # Test name 2
    assert config['groundstate/LiF'] == expected_config['groundstate/LiF']

    # Test name 3
    assert config['BSE/LiF'] == expected_config['BSE/LiF']

    return config


def test_configure_all_tests(defaults_config_dict, tests_config_dict):
    """
    Take a YAML-formatted string, parse and verify output dictionary
    (without removing any tests)

    Can't easily test `configure_all_tests` because it interacts with the file system
    As such, I've included the main calls here
    """
    test_names = [name for name in tests_config_dict.keys()]
    default_files_under_test = defaults_config_dict['default_files_under_test']

    # Result of routine that interacts with the file system
    default_inputs = {'groundstate/LDA_PW-PbTiO3':
                          ["INFO.OUT", "evalcore.xml", "geometry.xml", "eigval.xml", "atoms.xml"],
                      'groundstate/LiF':
                          ["INFO.OUT", "evalcore.xml", "geometry.xml", "eigval.xml", "atoms.xml"],
                      'BSE/LiF':
                          ['some_file']
                      }

    initialised_tests = initialise_tests(test_names, default_inputs, default_files_under_test)
    tests = assign_config_keys(initialised_tests, tests_config_dict)

    expected_tests = {'groundstate/LDA_PW-PbTiO3':
                          {'files_under_test': ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml'],
                           'inputs': ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml'],
                           'repeat': False,
                           'failing_builds': [CompilerBuild('intel', 'mpiandsmp'), CompilerBuild('intel', 'serial')],
                           'group': Group.SLOW_TESTS,
                           'comments': 'Some details on failing test to print from CI',
                           'depends_on': None,
                           'cmd_line_args': ''
                           },
                      'groundstate/LiF':
                          {'files_under_test': ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml'],
                           'inputs': ['INFO.OUT', 'evalcore.xml', 'geometry.xml', 'eigval.xml', 'atoms.xml'],
                           'repeat': False,
                           'failing_builds': [],
                           'group': Group.NONE,
                           'comments': '',
                           'depends_on': None,
                           'cmd_line_args': ''
                           },
                      'BSE/LiF':
                          {'files_under_test': ['some_file'],
                           'inputs': ['input.xml', 'STATE.OUT', 'Li.xml', 'F.xml'],
                           'repeat': True,
                           'failing_builds': [],
                           'group': Group.NONE,
                           'comments': '',
                           'depends_on': 'groundstate/LiF',
                           'cmd_line_args': ''
                           }
                      }

    # Compare values of CompilerBuild instances
    expected_PbTiO3_fbuilds = expected_tests['groundstate/LDA_PW-PbTiO3'].pop('failing_builds')
    PbTiO3_fbuilds = tests['groundstate/LDA_PW-PbTiO3'].pop('failing_builds')
    assert PbTiO3_fbuilds[0].compiler == expected_PbTiO3_fbuilds[0].compiler
    assert PbTiO3_fbuilds[1].compiler == expected_PbTiO3_fbuilds[1].compiler
    assert PbTiO3_fbuilds[0].build == expected_PbTiO3_fbuilds[0].build
    assert PbTiO3_fbuilds[1].build == expected_PbTiO3_fbuilds[1].build

    # Compare all other dict values
    assert tests['groundstate/LDA_PW-PbTiO3'] == expected_tests['groundstate/LDA_PW-PbTiO3']
    assert tests['groundstate/LiF'] == expected_tests['groundstate/LiF']
    assert tests['BSE/LiF'] == expected_tests['BSE/LiF']


def test_sort_tests_by_group():

    yaml_specified_tests = {
     'groundstate/LDA_PW-PbTiO3': {'group': Group.SLOW_TESTS},
     'groundstate/LiF': {'group': Group.NONE},
     'BSE/LiF': {'group': Group.NONE}
    }

    tests_by_group = sort_tests_by_group(yaml_specified_tests)

    expected_test_grouping = {enum: [] for enum in Group}
    expected_test_grouping.update({Group.NONE: ['groundstate/LiF', 'BSE/LiF'],
                                   Group.SLOW_TESTS: ['groundstate/LDA_PW-PbTiO3'],
                                   Group.SIRIUS: []
                                  })

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
