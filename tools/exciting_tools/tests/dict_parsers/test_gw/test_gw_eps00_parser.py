""" Test all GW output file parsers, except GW_INFO.OUT
"""
import pytest
import numpy as np

from excitingtools.utils.test_utils import MockFile
from excitingtools.utils.utils import get_new_line_indices

from excitingtools.exciting_dict_parsers.gw_eps00_parser import parse_eps00_frequencies, parse_eps00_gw, \
    _file_name


@pytest.fixture
def eps00_mock(tmp_path):
    # EPS00_GW.OUT segment. Note that first line of file must be blank.
    eps_string = """
(dielectric tensor, random phase approximation)

frequency index and value:      1    0.00529953
real part, imaginary part below
   8.31881773    0.00000000    0.00000000         0.00000000    0.00000000    0.00000000
   0.00000000    8.31881773    0.00000000         0.00000000    0.00000000    0.00000000
   0.00000000    0.00000000    8.31881773         0.00000000    0.00000000    0.00000000

frequency index and value:      2    0.02771249
real part, imaginary part below
   8.22189228    0.00000000    0.00000000        -0.00000000    0.00000000    0.00000000
   0.00000000    8.22189228    0.00000000         0.00000000   -0.00000000    0.00000000
   0.00000000    0.00000000    8.22189228         0.00000000    0.00000000   -0.00000000

frequency index and value:      3    0.06718440
real part, imaginary part below
   7.78004308    0.00000000    0.00000000         0.00000000    0.00000000    0.00000000
   0.00000000    7.78004308    0.00000000         0.00000000    0.00000000    0.00000000
   0.00000000    0.00000000    7.78004308         0.00000000    0.00000000    0.00000000

"""
    eps00_file = tmp_path / _file_name
    eps00_file.write_text(eps_string)
    return MockFile(eps00_file, eps_string)


def test_parse_eps00_frequencies(eps00_mock):
    """ Test parsing frequencies from EPS00_GW.OUT
    """
    line = get_new_line_indices(eps00_mock.string)

    ref = {1: 0.00529953, 2: 0.02771249, 3: 0.06718440}

    assert eps00_mock.string[line[0]:line[1]].isspace(), "First line of eps_string must be a whiteline"
    assert parse_eps00_frequencies(eps00_mock.string) == ref, "Frequency grid for eps00"


def test_parse_eps00_gw(eps00_mock):
    """ Test parsing EPS00_GW.OUT
    """
    line = get_new_line_indices(eps00_mock.string)
    assert eps00_mock.string[line[0]:line[1]].isspace(), "First line of eps_string must be a whiteline"

    ref = {
        1: {
            'frequency': 0.00529953,
            'eps00': {
                're': np.array([[8.31881773, 0., 0.], [0., 8.31881773, 0.], [0., 0., 8.31881773]]),
                'img': np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
                }
            },
        2: {
            'frequency': 0.02771249,
            'eps00': {
                're': np.array([[8.22189228, 0., 0.], [0., 8.22189228, 0.], [0., 0., 8.22189228]]),
                'img': np.array([[-0., 0., 0.], [0., -0., 0.], [0., 0., -0.]])
                }
            },
        3: {
            'frequency': 0.06718440,
            'eps00': {
                're': np.array([[7.78004308, 0.0, 0.0], [0., 7.78004308, 0.], [0., 0., 7.78004308]]),
                'img': np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
                }
            }
        }

    output = parse_eps00_gw(eps00_mock.file)

    assert len(output) == 3, "3 frequency points"
    assert [k for k in output[1].keys()] == ['frequency', 'eps00'], "Frequency point 1 keys "
    assert [k for k in output[2].keys()] == ['frequency', 'eps00'], "Frequency point 2 keys "
    assert [k for k in output[3].keys()] == ['frequency', 'eps00'], "Frequency point 3 keys "

    assert output[1]['frequency'] == 0.00529953, "Frequency point 1 value"
    assert output[2]['frequency'] == 0.02771249, "Frequency point 2 value"
    assert output[3]['frequency'] == 0.06718440, "Frequency point 3 value"

    assert np.allclose(output[1]['eps00']['re'], ref[1]['eps00']['re']),"Re{eps00} at frequency point 1"
    assert np.allclose(output[1]['eps00']['img'], ref[1]['eps00']['img']),"Im{eps00} at frequency point 1"

    assert np.allclose(output[2]['eps00']['re'], ref[2]['eps00']['re']),"Re{eps00} at frequency point 2"
    assert np.allclose(output[2]['eps00']['img'], ref[2]['eps00']['img']),"Im{eps00} at frequency point 2"

    assert np.allclose(output[3]['eps00']['re'], ref[3]['eps00']['re']),"Re{eps00} at frequency point 3"
    assert np.allclose(output[3]['eps00']['img'], ref[3]['eps00']['img']),"Im{eps00} at frequency point 3"
