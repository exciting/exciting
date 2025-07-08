import pytest
from excitingscripts.lattice.parameters import get_parameters
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def sgroup_out_mock(tmp_path) -> MockFile:
    """ Mock 'sgroup.out' file for testing get_parameters function. """
    sgroup_out_str = """Bravais lattice: Hexagonal

         a             b            c
     4.30000000    4.29999999   6.45000000
        alpha          beta          gamma
     90.00000000    90.00000000   120.00000011
    """

    sgroup_out_file = tmp_path / "sgroup.out"
    sgroup_out_file.write_text(sgroup_out_str)

    return MockFile(sgroup_out_file, sgroup_out_str)


def test_get_parameters(sgroup_out_mock, tmp_path):
    ref_results = {
        "sym": "Bravais lattice: Hexagonal",
        "a": 4.30000000,
        "b": 4.29999999,
        "c": 6.45000000,
        "alpha": 90.00000000,
        "beta": 90.00000000,
        "gamma": 120.00000011,
    }

    extracted_parameters = get_parameters(tmp_path)
    assert extracted_parameters == ref_results
