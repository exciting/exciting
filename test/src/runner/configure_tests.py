"""
Parse test suite configuration file, specified in YAML format, and use it in conjunction with default settings
to define the properties of all tests in the suite.
"""
from collections import namedtuple
from typing import List, Callable, Dict, Union
import enum
import yaml
import numpy as np
import os

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from .set_inputs import input_files_for_tests
from ..exciting_settings.constants import settings
from ..io.tolerances import list_tolerance_files_in_directory
from ..runner.profile import CompilerBuild
from ..tolerance.tol_classes import tol_file_name, tol_file_to_method
from ..tolerance.tol_classes import methods as valid_methods


class TestGroup(enum.Enum):
    """
    All GROUP definitions
    """
    NONE = enum.auto()
    SLOW_TESTS = enum.auto()


# Properties/attributes of a test case.
# TODO(Alex/Ben/Hannah) # Issue 114. Add 'depends_on'
ConfigurationDefaults = namedtuple('ConfigurationDefaults',
                                   ['files_under_test', 'inputs', 'repeat', 'failing_builds', 'group', 'comments']
                                   )


def str_attributes_to_enums(data: dict):
    """
    For test properties/attributes defined in input as strings, find and replace them with enums.

    For example, an input:
      {'groundstate/LDA_PW-PbTiO3': {'group': SLOW_TESTS}}

    would be returned as:

     {'groundstate/LDA_PW-PbTiO3': {'group': TestGroup.SLOW_TESTS}}

    :param dict data: Tests parsed from config, with values for `failing_builds` and `group` defined as strings.
    :return dict data: Tests parsed from config, with string values replaced with enums.
    """
    # For isinstance(value, str)
    def replace_scalar(key: str, str_to_enum: Callable) -> dict:
        for test_name in data.keys():
            try:
                str_value = data[test_name][key]
                data[test_name][key] = str_to_enum(str_value)
            except (KeyError, TypeError):
                # Only replace when the data is explicitly specified in the parsed config dict
                pass
        return data

    # For type(value) == List[str] or List[List[str]]
    def replace_list(key: str, str_to_enum: Callable) -> dict:
        for test_name in data.keys():
            try:
                list_value = data[test_name][key]
                data[test_name][key] = [str_to_enum(string) for string in list_value]
            except (KeyError, TypeError):
                # Only replace when the data is explicitly specified in the parsed config dict
                pass
        return data

    def compiler_build_str_to_enum(string: str):
        try:
            compiler, build = string.split('_')[0:2]
        except ValueError:
            raise ValueError(f'String representation of CompilerBuild enum is invalid: {string}. \n'
                             f'See src/runner/profile.py or config file documentation for valid strings.')

        return CompilerBuild(compiler, build)

    # `data` mutable and passed implicitly to these local routines, as all defined in same scope
    data = replace_scalar('group', lambda string: TestGroup[string])
    data = replace_list('failing_builds', compiler_build_str_to_enum)

    return data


def check_default_files_under_test(default_files_under_test: dict):
    """
    Check that all methods have "default files under test" specified in config.

    Expect an input of the form:
      default_files_under_test = {'groundstate': ["INFO.OUT","evalcore.xml", "geometry.xml", "eigval.xml", "atoms.xml"],
                                   ...
                                  'method': files_under_test_for_method
                                  }

    :param dict default_files_under_test: Files under test for all methods
    """
    tabulated_methods = {m.lower() for m in tol_file_name.keys()}
    methods_from_config = {m.lower() for m in default_files_under_test.keys()}

    missing_methods = tabulated_methods - methods_from_config
    if missing_methods:
        raise ValueError(f"Missing default_files_under_test in config, for methods: {missing_methods}")

    missing_methods = methods_from_config - tabulated_methods
    if missing_methods:
        raise ValueError(f"Missing tabulated methods in tol_classes.py, for methods: {missing_methods}")


def check_groups_in_config(groups: dict):
    """
    Check that all groups explicitly defined in the defaults config file are also defined as enums.

    If one is missing in either case, add it.

    Expect an input of the form:
        groups = {NONE: True, SLOW_TESTS: False, EXPERIMENTAL_TESTS: False, ... GROUP: to_run}

    :param dict groups: Groups read from config, with keys = group names and values = True/False.
    """
    enum_group_names = set(TestGroup._member_names_)
    groups_in_config_file = set(groups.keys())
    if enum_group_names != groups_in_config_file:
        raise ValueError('Groups specified in TestGroup and config file differ')


def parse_config_defaults_file(yaml_str: str) -> dict:
    """
    Read the defaults configuration file and extract data.

    :param str yaml_str: String with YAML formatting.
    :return dict config_defaults: Dict of defaults. Refer to the input file for the structure of each value.
    """
    config_defaults: dict = yaml.load(yaml_str, Loader=Loader)
    check_default_files_under_test(config_defaults['default_files_under_test'])
    check_groups_in_config(config_defaults['group_execution'])
    return config_defaults


def get_method(test_name: str) -> str:
    """
    Get method from the test name.

    TODO Remove interaction with the file system.
    One needs to know the method to choose the correct default `files_under_test`,
    when they are not specified in the config file.

    The difficultly is that in most test cases, the file path defines the method.
    For example, `BSE/LiF`. However, properties are all grouped under one subdirectory of the same name.
    In this case, one can infer from the method from the tolerance file.
    This is annoying to unit test though, as it requires one to interact with the file system. So every routine
    that uses `get_method` (i.e. all higher-level routines in configure_tests.py) will then require
    some file mocking when unit-tested.

    Hence this use of this solution of a) checking the test name string and b) falling back to the tolerance file iff
    it's a `properties` test.

    The other solutions one could implement are:
     a) Don't interact with the file system and instead specify every method in config input
        - This is the cleanest solution with minimal extra typing
        - May also restrict developers testing multiple properties in one test.
     b) Split up the properties subdirectory according to the different tolerance files.

    :param str test_name: Test name, prepended by path.
    :return str method: Method name.
    """
    method = get_method_from_test_name(test_name)
    if isinstance(method, KeyError):
        method = get_method_from_tolerance_file(test_name)
    return method


def get_method_from_test_name(test_name: str) -> Union[str, KeyError]:
    """
    Get method string from test name.

    Specific to the config file test naming convention.
    Expect test_name = 'method/test_case' or 'test_farm/method/test_case'
    Note, relying on strings and dealing with case is not the most robust approach.

    :param str test_name: Test name, prepended by path.
    :return Union[str, KeyError] method: Method name or KeyError
    """
    method = test_name.split('/')[-2].lower()
    if method not in valid_methods:
        # Properties are not differentiated in the directory structure of the test_farm, hence this requirement
        return KeyError()
    return method


def get_method_from_tolerance_file(test_name: str, subdirectory=settings.ref_dir) -> str:
    """
    Get method from tolerance file inspection.

    Uses the method specified in the {tolerance file : method map}

    :param str test_name: Test name, prepended by path.
    :param optional, str subdirectory: Subdirectory of test case, in which to look for tolerance file.
    :return str method: Method name.
    """
    directory = os.path.join(test_name, subdirectory)
    tolerance_files = list_tolerance_files_in_directory(directory)

    if len(tolerance_files) != 1:
        raise ValueError(f"Found {len(tolerance_files)} tolerance files in directory: {directory}")

    tol_file = tolerance_files[0]
    try:
        method = tol_file_to_method[tol_file]
    except KeyError:
        method_and_testname = "/".join(s for s in test_name.split('/')[-2:])
        valid_tol_files_str = "\n".join(m for m in tol_file_to_method.keys())
        raise KeyError(f'Invalid tolerance file, {tol_file}, for {method_and_testname}. \n'
                       f'Method must correspond to a tolerance:\n{valid_tol_files_str}'
                       )

    return method


def initialise_tests(test_names: List[str], default: ConfigurationDefaults) -> dict:
    """
    Create tests dictionary, with all tests initialised using default properties.

    :param List[str] test_names: List of test names.
    :param ConfigurationDefaults default: Default test attributes.
    :return dict tests: Initialised dictionary of test cases.
    """
    # defaults.inputs *currently* defined from searching test_farm.
    # Expect default.inputs = {'PBE-Al': ['input.xml', 'Ar.xml'], ..., 'test_case': ['input', 'files']}
    inconsistent_test_names = set(default.inputs.keys()) - set(test_names)
    if inconsistent_test_names:
        raise KeyError(f"Default input keys (from inspecting test_farm) inconsistent with test_names (from config file): "
                       f"{inconsistent_test_names}")

    # TODO(Alex/Ben/Hannah) # Issue 114. Add 'depends_on'
    tests = {}
    for name in test_names:
        method = get_method(name)
        default_files_under_test = default.files_under_test[method]
        tests[name] = {'files_under_test': default_files_under_test,
                       'inputs': default.inputs[name],
                       'repeat': default.repeat,
                       'failing_builds': default.failing_builds,
                        #'depends_on': default.depends_on,
                       'group': default.group,
                       'comments': default.comments
                       }

    # Check all keys are assigned default values.
    random_test_defaults = next(iter(tests.values()))
    assert set(random_test_defaults.keys()) == set(ConfigurationDefaults._fields), \
        "Field present in namedtuple is missing from this initialisation."

    return tests


def setup_tests(config_data: dict, default_attributes: ConfigurationDefaults) -> dict:
    """
    Initialise a dictionary of test cases with default properties, then replace the defaults
    with explicitly-specified attributes, when specified in the tests config file.

    :param dict config_data: Tests from the config file, with explicitly-defined properties.
    :param ConfigurationDefaults default_attributes: Default properties/attributes.
    :return dict tests: All tests in the suite, with defaults overwritten with config file properties,
    where explicitly defined.
    """
    test_names = [name for name in config_data.keys()]
    initialised_tests = initialise_tests(test_names, default_attributes)
    tests = assign_test_configurations(initialised_tests, config_data)
    return tests


def check_tests_have_all_configuration_fields(tests: dict):
    """
    Check that each test case has all properties defined by the Configuration tuple.

    :param dict tests: Test dictionary, with initialised values.
    """
    expected_keys = set(ConfigurationDefaults._fields)

    for name, attributes in tests.items():
        test_keys = set(attributes.keys())
        if expected_keys != test_keys:
            raise KeyError(f'Test dictionary {name}, not initialised with all ConfigurationDefaults keys')


def assign_test_configurations(initialised_tests: dict, test_configs: dict) -> dict:
    """
    Add any attributes from tests with explicitly-defined attributes defined in the config file,
    overwriting the defaults.

    initialised_tests gets mutated.

    :param dict initialised_tests: Test dictionary, with initialised values.
    :param dict test_configs: Test data parsed from config file
    :return dict tests: Test dictionary, with any attributes defined in test_configs used to overwrite the defaults.
    """
    # Should be guaranteed by `initialise_tests`, hence maybe remove
    check_tests_have_all_configuration_fields(initialised_tests)

    for name, attributes in test_configs.items():
        try:
            for key, value in attributes.items():
                initialised_tests[name][key] = value
        except AttributeError:
            pass
    return initialised_tests


def configure_all_tests(tests_yaml: str, defaults_yaml: str) -> dict:
    """
    Configure all tests in the suite.

    Parse two config files containing default settings, and test configurations, initialise tests with default
    settings, then replace any properties that have been explicitly specified in the test config file.

    :param str tests_yaml: Test config file as a string.
    :param str defaults_yaml: Defaults config file as a string.
    :return dict tests: Tests dict, containing all test cases (even those that shouldn't run) and their properties.
    """
    config_defaults = parse_config_defaults_file(defaults_yaml)

    config_data: dict = yaml.load(tests_yaml, Loader=Loader)
    config_data = str_attributes_to_enums(config_data)
    # Must prepend with 'test_farm' directory, as this is not specified in the naming convention of the config file
    config_data = {os.path.join(settings.test_farm, name): value for name, value in config_data.items()}

    test_names = [name for name in config_data.keys()]
    default_inputs = input_files_for_tests(test_names, subdirectory='ref')
    # TODO(Alex/Ben/Hannah) # Issue 114. Add 'depends_on'
    default_attributes = ConfigurationDefaults(files_under_test=config_defaults['default_files_under_test'],
                                               inputs=default_inputs,
                                               repeat=False,
                                               failing_builds=[],
                                               # depends_on=[],
                                               group=TestGroup.NONE,
                                               comments=''
                                               )

    tests = setup_tests(config_data, default_attributes)
    return tests


def find_tests_with_attribute(tests: dict, attribute_key: str, default_value) -> dict:
    """
    Create a dictionary of tests with a specific attribute, as specified by attribute_key.

    For example, attribute_key = failing tests
    return failing_tests = {'test1': failing_builds_i, 'testi': failing_builds_1}

    :param dict tests: Tests.
    :param str attribute_key: Key of attribute.
    :param default_value: Default value of attribute, as specified in 'default'.
    :return dict: subset of tests.
    """
    test_names = []
    specific_attributes = []

    for name, attributes in tests.items():
        test_names.append(name)
        specific_attributes.append(attributes[attribute_key])

    not_equal = [entry != default_value for entry in specific_attributes]
    indices = np.where(not_equal)[0]
    test_names_with_attribute = [test_names[i] for i in indices]

    return {name: tests[name] for name in test_names_with_attribute}


def sort_tests_by_group(tests: dict) -> Dict[TestGroup, List[str]]:
    """
    Sort tests according to their `group` tag/attribute.

    For example:

    {SLOW_TESTS: ['test_farm/XANES/TiO2'],
     GROUP_NAME: ['test_farm/method/test_name', 'test_farm/method/another_test']
    }

    :param dict tests: Tests.
    :return Dict[TestGroup, List[str]] tests_by_group: Tests sorted by group.
    """
    # Initialise
    tests_by_group = {group: [] for group in TestGroup}

    # Log tests by 'group' attribute
    for name, attributes in tests.items():
        group = attributes['group']
        tests_by_group[group].append(name)

    return tests_by_group


def find_tests_to_skip_by_group(all_tests: dict, group_execution: dict) -> set:
    """
    List tests with groups that are not assigned to run.

    :param dict all_tests: All test cases.
    :param dict group_execution: All groups, and whether they should run.
     For example {SLOW_TESTS: False, GROUP_NAME: True, ...}

    :return dict tests_to_skip: Tests to skip.
    """
    tests_by_group = sort_tests_by_group(all_tests)

    tests_to_run = set()
    for group_str, to_run in group_execution.items():
        if to_run:
            group = TestGroup[group_str]
            tests_to_run.update(tests_by_group[group])

    tests_to_skip = set(all_tests.keys()) - tests_to_run
    return tests_to_skip
