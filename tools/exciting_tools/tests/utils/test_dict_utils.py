"""
Tests for functions in dict_utils.py
"""
import numpy as np

from excitingtools.utils.dict_utils import container_converter, serialise_dict_values, delete_nested_key


def test_convert_container():
    """Test container_converter function on string values in a dict."""
    input = {
        'a': '5.0',
        'b': ['1.0', '10.0'],
        'c': {
            'a1': '1.0'
        },
        'd': "blu",
        'e': [['1', '2.3'], '3'],
        'f': ['1', 'a']
    }
    expected = {
        'a': 5.0,
        'b': [1.0, 10.0],
        'c': {
            'a1': 1.0
        },
        'd': "blu",
        'e': [[1, 2.3], 3],
        'f': [1, 'a']
    }

    assert container_converter(input) == expected, (
        'String value/s failed to convert to numerical values')


def test_convert_container_no_strings():
    """Test container_converter where there are no string values in dict."""
    input = {
        'a': 5.0,
        'b': [1.0, 10.0],
        'c': {
            'a1': 1.0
        },
        'd': "blu",
        'e': [[1, 2.3], 3],
        'f': [1, 11.3]
    }

    assert container_converter(input) == input, (
        "Expect the converter to do nothing")


def test_convert_container_data_type():
    input = {
        'a': '5.0',
        'b': ['1.0', '10.0'],
        'c': {
            'a1': '1.0'
        },
        'd': "blu",
        'e': [['1', '2.3'], '3']
    }
    expected = {
        'a': 5.0,
        'b': [1.0, 10.0],
        'c': {
            'a1': 1.0
        },
        'd': "blu",
        'e': [[1, 2.3], 3]
    }

    output = container_converter(input)

    for elem1, elem2 in zip(output.values(), expected.values()):
        assert type(elem1) == type(elem2)
    for elem, elem2 in zip(output['c'].values(), expected['c'].values()):
        assert type(elem) == type(elem2)
    for elem, elem2 in zip(output['b'], expected['b']):
        assert type(elem) == type(elem2)
    for elem, elem2 in zip(output['e'], expected['e']):
        assert type(elem) == type(elem2)
    for elem, elem2 in zip(output['e'][0], expected['e'][0]):
        assert type(elem) == type(elem2)

    assert output == expected, (
        "Output is consistent with the expected dictionary")


class Mock:

    def __init__(self, a, b):
        self.a = a
        self.b = b


def test_serialise_dict_values():

    # Value is an object
    input = {'mock_key': Mock(1, 2)}
    output = serialise_dict_values(input)
    assert output == {
        'mock_key': {
            'a': 1,
            'b': 2
        }
    }, "Convert object values to dicts"

    # Object nested in a dictionary
    input = {'mock_key': {'another-key': Mock(1, 2)}}
    output = serialise_dict_values(input)
    assert output == {
        'mock_key': {
            'another-key': {
                'a': 1,
                'b': 2
            }
        }
    }, "Convert nested object values into dicts"

    # Object nested in a list, and within a dictionary within a list
    input = {'a': [1, 2, Mock(3, 4), {'b': Mock(5, 6)}]}
    output = serialise_dict_values(input)
    assert output == {'a': [1, 2, {'a': 3, 'b': 4}, {'b': {'a': 5, 'b': 6}}]}, \
        "Convert nested object values into dicts, where the top-level container value is a list"


def test_serialise_dict_value_is_tuple():

    # Object nested in a dictionary, nested within a tuple
    input = {'mock_key': (1, {'another-key': Mock(1, 2)})}
    output = serialise_dict_values(input)
    assert output == {'mock_key': [1, {'another-key': {'a': 1, 'b': 2}}]}, \
        "Convert nested object values into dicts, where the top-level container value is a tuple." \
        "Note, tuple cannot be mutated so it is converted to a list by serialise_dict_values"


def test_serialise_dict_value_is_set():

    # Object nested in a dictionary, nested within a set
    input = {'mock_key': {1, 2, Mock(1, 2)}}
    output = serialise_dict_values(input)

    output_keys = [k for k in output.keys()]
    output_values = [v for v in output.values()]

    assert len(output_keys) == 1
    assert output_keys[0] == 'mock_key'
    assert len(output_values) == 1

    # Mock object -> dictionary, meaning one cannot compare the input and outout with collections.Counter,
    # sets or sorted (neither object nor dictionary is hashable)
    def unordered_lists_the_same(a: list, b: list) -> bool:
        """
        Check contents of two lists are consistent, irregardless of order
        """
        for element in a:
            try:
                b.remove(element)
            except ValueError:
                return False
        return not b

    assert unordered_lists_the_same(output_values[0], [1, 2, {'a': 1, 'b': 2}]), \
        "Convert nested object values into dicts, where the top-level container value is a set. " \
        "Note, set is not iterable and is converted to a list with some arbitrary order by serialise_dict_values"


def test_serialise_dict_values_null_behaviour():
    """"
    Test serialize_dict_values on different dict key types.
    """
    input = {'mock_key': [1, 2]}
    output = serialise_dict_values(input)
    assert output == input, "Pass over list values"

    input = {'mock_key': (1, 2)}
    output = serialise_dict_values(input)
    assert output == input, "Pass over tuple values"

    input = {'mock_key': {1, 2}}
    output = serialise_dict_values(input)
    assert output == input, "Pass over set values"

    input = {'mock_key': np.array([1., 2., 3.])}
    output = serialise_dict_values(input)
    assert output == input, "Pass over np.array values"


def test_delete_nested_key():
    """
    Test delete_nested_key removes a key from dict.
    """
    input = {'a': {'b': 1, 'c': {'d': 1, 'e': 2}}}
    key_to_remove = ['a']
    delete_nested_key(input, key_to_remove)
    output = {}
    assert output == input

    input = {'a': {'b': 1, 'c': {'d': 1, 'e': 2}}}
    key_to_remove = ['a', 'c', 'd']
    delete_nested_key(input, key_to_remove)
    output = {'a': {'b': 1, 'c': {'e': 2}}}
    assert output == input
