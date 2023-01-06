import numpy as np
import pytest
from excitingtools.utils.test_utils import MockFile
from excitingtools.exciting_obj_parsers.eigenvalue_parser import parse_eigenvalues


@pytest.fixture
def eigval_xml_mock(tmp_path) -> MockFile:
    """ Mock 'eigval.xml' data, containing only two k-sampling points.
    """
    eigval_xml_str = """<?xml version="1.0" encoding="UTF-8"?>
    <eigval>
        <kpt ik="1" vkl="0.000000000000e0 0.000000000000e0 0.000000000000e0">
            <state ist="1" eigenvalue="-3.728453763676e-1" occupancy="2.000000000000e0"/>
            <state ist="2" eigenvalue="4.171109202867e-1" occupancy="2.000000000000e0"/>
            <state ist="3" eigenvalue="4.171125401350e-1" occupancy="2.000000000000e0"/>
            <state ist="4" eigenvalue="4.171125401350e-1" occupancy="2.000000000000e0"/>
            <state ist="5" eigenvalue="6.240763860566e-1" occupancy="0.000000000000e0"/>
            <state ist="6" eigenvalue="6.240773264914e-1" occupancy="0.000000000000e0"/>
            <state ist="7" eigenvalue="6.240773264914e-1" occupancy="0.000000000000e0"/>
            <state ist="8" eigenvalue="9.057770677021e-1" occupancy="0.000000000000e0"/>
            <state ist="9" eigenvalue="1.129829432130e0" occupancy="0.000000000000e0"/>
            <state ist="10" eigenvalue="1.392608317584e0" occupancy="0.000000000000e0"/>
        </kpt>
        <kpt ik="2" vkl="2.500000000000e-1 0.000000000000e0 0.000000000000e0">
            <state ist="1" eigenvalue="-3.098896419810e-1" occupancy="2.000000000000e0"/>
            <state ist="2" eigenvalue="1.757980236173e-1" occupancy="2.000000000000e0"/>
            <state ist="3" eigenvalue="3.536527660909e-1" occupancy="2.000000000000e0"/>
            <state ist="4" eigenvalue="3.536539162999e-1" occupancy="2.000000000000e0"/>
            <state ist="5" eigenvalue="7.055469643798e-1" occupancy="0.000000000000e0"/>
            <state ist="6" eigenvalue="7.055481190802e-1" occupancy="0.000000000000e0"/>
            <state ist="7" eigenvalue="7.176489731254e-1" occupancy="0.000000000000e0"/>
            <state ist="8" eigenvalue="9.578009362424e-1" occupancy="0.000000000000e0"/>
            <state ist="9" eigenvalue="1.183488522547e0" occupancy="0.000000000000e0"/>
            <state ist="10" eigenvalue="1.242232237157e0" occupancy="0.000000000000e0"/>
        </kpt>
    </eigval> 
    """
    eigval_xml_file = tmp_path / "eigval.xml"
    eigval_xml_file.write_text(eigval_xml_str)
    return MockFile(eigval_xml_file, eigval_xml_str)


def test_parse_eigenvalues(eigval_xml_mock):
    eigval_data = parse_eigenvalues(eigval_xml_mock.file)

    ref_k_points = np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]])

    ref_eigenvalues = np.array([[-0.37284538, 0.41711092, 0.41711254, 0.41711254, 0.62407639, 0.62407733, 0.62407733,
                                 0.90577707, 1.12982943, 1.39260832],
                                [-0.30988964, 0.17579802, 0.35365277, 0.35365392, 0.70554696, 0.70554812, 0.71764897,
                                 0.95780094, 1.18348852, 1.24223224]])

    ref_occupations = np.array([[2, 2, 2, 2, 0, 0, 0, 0, 0, 0], [2, 2, 2, 2, 0, 0, 0, 0, 0, 0]])

    assert eigval_data.state_range.first_state == 1
    assert eigval_data.state_range.last_state == 10
    assert np.allclose(eigval_data.k_points, ref_k_points)
    assert eigval_data.k_indices == [1, 2]
    assert np.allclose(eigval_data.all_eigenvalues, ref_eigenvalues)
    assert np.allclose(eigval_data.occupations, ref_occupations)
