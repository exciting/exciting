import pytest
import numpy as np

from excitingtools.utils.test_utils import MockFile
from excitingtools.exciting_dict_parsers.gw_vxc_parser import parse_vxcnn, vkl_from_vxc


@pytest.fixture
def vxc_mock(tmp_path):
    """ Mock VXCNN.DAT data with energies for 3 k-points
    """
    vxc_string = """ik=   1    vkl=  0.0000  0.0000  0.0000
       1       -2.908349       -0.000000
       2       -2.864103        0.000000
       3       -2.864103       -0.000000
       4       -2.864103       -0.000000
       5       -2.764246       -0.000000
       6       -2.764246        0.000000
       7       -2.764809       -0.000000
       8       -2.764809       -0.000000
       9       -2.764809       -0.000000
      10       -1.034312       -0.000000
      11       -0.929773       -0.000000

    ik=   2    vkl=  0.0000  0.0000  0.5000
     1       -2.908349        0.000000
     2       -2.864100       -0.000000
     3       -2.864100        0.000000
     4       -2.864101       -0.000000
     5       -2.764227       -0.000000
     6       -2.764227        0.000000
     7       -2.764770        0.000000
     8       -2.764777        0.000000
     9       -2.764777        0.000000
    10       -1.037228       -0.000000
    11       -0.889360        0.000000

    ik=   3    vkl=  0.0000  0.5000  0.5000
       1       -2.908349       -0.000000
       2       -2.864099       -0.000000
       3       -2.864099       -0.000000
       4       -2.864099        0.000000
       5       -2.764195        0.000000
       6       -2.764208        0.000000
       7       -2.764760       -0.000000
       8       -2.764780       -0.000000
       9       -2.764780       -0.000000
      10       -1.038185        0.000000
      11       -0.887485       -0.000000
    """
    file = tmp_path / "VXCNN.DAT"
    file.write_text(vxc_string)
    return MockFile(file, vxc_string)


def test_vkl_from_vxc(vxc_mock):
    # k-points in units of the lattice vectors
    vkl_ref = {
        1: [0.0000, 0.0000, 0.0000], 2: [0.0000, 0.0000, 0.5000], 3: [0.0000, 0.5000, 0.5000]
        }
    output = vkl_from_vxc(vxc_mock.string)
    assert len(output) == 3, "Expect 3 k-points"
    assert output == vkl_ref, "vkl values equal to vkl_ref"


def test_parse_vxcnn(vxc_mock):

    # Reference V_xc extracted from vxc_string, defined above
    v_xc_1 = np.array([[-2.908349, -0.000000], [-2.864103, 0.000000], [-2.864103, -0.000000],
                       [-2.864103, -0.000000], [-2.764246, -0.000000], [-2.764246, 0.000000],
                       [-2.764809, -0.000000], [-2.764809, -0.000000], [-2.764809, -0.000000],
                       [-1.034312, -0.000000], [-0.929773, -0.000000]])

    v_xc_2 = np.array([[-2.908349, 0.000000], [-2.864100, -0.000000], [-2.864100, 0.000000],
                       [-2.864101, -0.000000], [-2.764227, -0.000000], [-2.764227, 0.000000],
                       [-2.764770, 0.000000], [-2.764777, 0.000000], [-2.764777, 0.000000],
                       [-1.037228, -0.000000], [-0.889360, 0.000000]])

    v_xc_3 = np.array([[-2.908349, -0.000000], [-2.864099, -0.000000], [-2.864099, -0.000000],
                       [-2.864099, 0.000000], [-2.764195, 0.000000], [-2.764208, 0.000000],
                       [-2.764760, -0.000000], [-2.764780, -0.000000], [-2.764780, -0.000000],
                       [-1.038185, 0.000000], [-0.887485, -0.000000]])

    output = parse_vxcnn(vxc_mock.file)

    assert [key for key in output[1].keys()] == ['vkl', 'v_xc_nn'], "Key consistency for ik=1 of parsed vxcnn"
    assert [key for key in output[2].keys()] == ['vkl', 'v_xc_nn'], "Key consistency for ik=2 of parsed vxcnn"
    assert [key for key in output[3].keys()] == ['vkl', 'v_xc_nn'], "Key consistency for ik=3 of parsed vxcnn"

    assert output[1]['vkl'] == [0.0000, 0.0000, 0.0000], "vkl (ik=1)"
    assert output[2]['vkl'] == [0.0000, 0.0000, 0.5000], "vkl (ik=2)"
    assert output[3]['vkl'] == [0.0000, 0.5000, 0.5000], "vkl (ik=3)"

    assert output[1]['v_xc_nn'].shape == (11, 2), "Expect V_xc to have 2 cols for 11 states"
    assert output[2]['v_xc_nn'].shape == (11, 2), "Expect V_xc to have 2 cols for 11 states"
    assert output[3]['v_xc_nn'].shape == (11, 2), "Expect V_xc to have 2 cols for 11 states"

    assert np.allclose(output[1]['v_xc_nn'], v_xc_1), "v_xc_nn for ik=1"
    assert np.allclose(output[2]['v_xc_nn'], v_xc_2), "v_xc_nn for ik=2"
    assert np.allclose(output[3]['v_xc_nn'], v_xc_3), "v_xc_nn for ik=3"
