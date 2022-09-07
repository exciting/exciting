import numpy as np
import pytest
from excitingtools.utils.test_utils import MockFile
from excitingtools.exciting_dict_parsers.properties_parser import parse_band_structure_xml, parse_band_structure_dat


@pytest.fixture
def band_structure_dat_mock(tmp_path) -> MockFile:
    """ Mock 'bandstructure.dat' data, containing only two bands and
    only 6 k-sampling points per band.
    """
    bs_dat_str = """#            1          2         6
     1     1    1.000000    0.000000    0.000000   0.000000000      -3.370713328    
     1     2    0.988281    0.011719    0.000000  0.2697635234E-01  -3.370710744    
     1     3    0.976562    0.023438    0.000000  0.5395270469E-01  -3.370705193    
     1     4    0.964844    0.035156    0.000000  0.8092905703E-01  -3.370698602    
     1     5    0.953125    0.046875    0.000000  0.1079054094      -3.370682200    
     1     6    0.941406    0.058594    0.000000  0.1348817617      -3.370661229   

     2     1    1.000000    0.000000    0.000000   0.000000000      -2.024168147    
     2     2    0.988281    0.011719    0.000000  0.2697635234E-01  -2.024186985    
     2     3    0.976562    0.023438    0.000000  0.5395270469E-01  -2.024297489    
     2     4    0.964844    0.035156    0.000000  0.8092905703E-01  -2.024460642    
     2     5    0.953125    0.046875    0.000000  0.1079054094      -2.024597185    
     2     6    0.941406    0.058594    0.000000  0.1348817617      -2.024765908    
    """
    bs_dat_file = tmp_path / "bandstructure.dat"
    bs_dat_file.write_text(bs_dat_str)
    return MockFile(bs_dat_file, bs_dat_str)


def test_parse_band_structure_xml():
    band_data = parse_band_structure_xml(band_structure_xml)

    assert band_data['n_kpts'] == band_data['band_energies'].shape[0], (
        "First dim of bands array equals the number of k-sampling points in the band structure")
    assert band_data['n_kpts'] == 6, "sampling points per band"
    assert band_data['n_bands'] == 2, "band_structure_xml contains two bands"

    ref_k_points = [0., 0.04082159, 0.08164318, 0.12246477, 0.16328636, 0.20410795]

    ref_bands = np.array([[-0.45003454, -0.00937631],
                          [-0.44931675, -0.01419609],
                          [-0.44716535, -0.02681183],
                          [-0.44358550, -0.04401707],
                          [-0.43858583, -0.06361773],
                          [-0.43217908, -0.0844593]])

    assert np.allclose(band_data['k_points_along_band'], ref_k_points, atol=1.e-8)
    assert np.allclose(band_data['band_energies'], ref_bands, atol=1.e-8)


def test_parse_band_structure_xml_vertices():
    vertices_ref = [{'distance': 0.0, 'label': 'G', 'coord': [0.0, 0.0, 0.0]},
                    {'distance': 0.6123238446, 'label': 'X', 'coord': [0.5, 0.0, 0.5]},
                    {'distance': 0.918485767, 'label': 'W', 'coord': [0.5, 0.25, 0.75]},
                    {'distance': 1.134974938, 'label': 'K', 'coord': [0.375, 0.375, 0.75]},
                    {'distance': 1.784442453, 'label': 'G', 'coord': [0.0, 0.0, 0.0]},
                    {'distance': 2.314730457, 'label': 'L', 'coord': [0.5, 0.5, 0.5]},
                    {'distance': 2.689700702, 'label': 'U', 'coord': [0.625, 0.25, 0.625]},
                    {'distance': 2.906189873, 'label': 'W', 'coord': [0.5, 0.25, 0.75]},
                    {'distance': 3.339168216, 'label': 'L', 'coord': [0.5, 0.5, 0.5]},
                    {'distance': 3.71413846, 'label': 'K', 'coord': [0.375, 0.375, 0.75]},
                    {'distance': 3.71413846, 'label': 'U', 'coord': [0.625, 0.25, 0.625]},
                    {'distance': 3.930627631, 'label': 'X', 'coord': [0.5, 0.0, 0.5]}]

    band_data = parse_band_structure_xml(band_structure_xml)
    assert band_data['vertices'] == vertices_ref


def test_parse_band_structure_dat(band_structure_dat_mock):
    band_data = parse_band_structure_dat(band_structure_dat_mock.file)

    assert band_data['n_kpts'] == band_data['band_energies'].shape[0], (
        "First dim of bands array equals the number of k-sampling points in the band structure")
    assert band_data['n_kpts'] == 6, "sampling points per band"
    assert band_data['n_bands'] == 2, "band_structure_xml contains two bands"

    ref_k_points = np.array([[1.000000, 0.000000, 0.000000],
                             [0.988281, 0.011719, 0.000000],
                             [0.976562, 0.023438, 0.000000],
                             [0.964844, 0.035156, 0.000000],
                             [0.953125, 0.046875, 0.000000],
                             [0.941406, 0.058594, 0.000000]])

    ref_bands = np.array([[-3.370713328, -2.024168147],
                          [-3.370710744, -2.024186985],
                          [-3.370705193, -2.024297489],
                          [-3.370698602, -2.024460642],
                          [-3.370682200, -2.024597185],
                          [-3.370661229, -2.024765908]])

    ref_flattened_k_points = np.array([0., 0.02697635, 0.05395270, 0.08092905, 0.10790540, 0.13488176])

    assert np.allclose(band_data['k_points'], ref_k_points, atol=1.e-8)
    assert np.allclose(band_data['flattened_k_points'], ref_flattened_k_points, atol=1.e-8)
    assert np.allclose(band_data['band_energies'], ref_bands, atol=1.e-8)

# Band structure of silicon, containing two lowest bands and only 6 k-sampling points per band

band_structure_xml = """<?xml version="1.0" encoding="UTF-8"?>
<?xml-stylesheet href="http://xml.exciting-code.org/visualizationtemplates/bandstructure2html.xsl" type="text/xsl"?>
<bandstructure>
  <title>silicon-primitive-PBEsol</title>
  <band>
    <point distance="0.000000000" eval="-0.4500345371"/>
    <point distance="0.4082158964E-01" eval="-0.4493167456"/>
    <point distance="0.8164317929E-01" eval="-0.4471653536"/>
    <point distance="0.1224647689" eval="-0.4435855004"/>
    <point distance="0.1632863586" eval="-0.4385858287"/>
    <point distance="0.2041079482" eval="-0.4321790844"/>
  </band>
  <band>
    <point distance="0.000000000" eval="-0.9376311739E-02"/>
    <point distance="0.4082158964E-01" eval="-0.1419608531E-01"/>
    <point distance="0.8164317929E-01" eval="-0.2681183031E-01"/>
    <point distance="0.1224647689" eval="-0.4401706572E-01"/>
    <point distance="0.1632863586" eval="-0.6361773499E-01"/>
    <point distance="0.2041079482" eval="-0.8445929704E-01"/>
  </band>
  <vertex distance="0.000000000" upperboundary="0.8623187822" lowerboundary="-1.106211197" label="G" coord="0.000000000       0.000000000       0.000000000"/>
  <vertex distance="0.6123238446" upperboundary="0.8623187822" lowerboundary="-1.106211197" label="X" coord="0.5000000000       0.000000000      0.5000000000"/>
  <vertex distance="0.9184857670" upperboundary="0.8623187822" lowerboundary="-1.106211197" label="W" coord="0.5000000000      0.2500000000      0.7500000000"/>
  <vertex distance="1.134974938" upperboundary="0.8623187822" lowerboundary="-1.106211197" label="K" coord="0.3750000000      0.3750000000      0.7500000000"/>
  <vertex distance="1.784442453" upperboundary="0.8623187822" lowerboundary="-1.106211197" label="G" coord="0.000000000       0.000000000       0.000000000"/>
  <vertex distance="2.314730457" upperboundary="0.8623187822" lowerboundary="-1.106211197" label="L" coord="0.5000000000      0.5000000000      0.5000000000"/>
  <vertex distance="2.689700702" upperboundary="0.8623187822" lowerboundary="-1.106211197" label="U" coord="0.6250000000      0.2500000000      0.6250000000"/>
  <vertex distance="2.906189873" upperboundary="0.8623187822" lowerboundary="-1.106211197" label="W" coord="0.5000000000      0.2500000000      0.7500000000"/>
  <vertex distance="3.339168216" upperboundary="0.8623187822" lowerboundary="-1.106211197" label="L" coord="0.5000000000      0.5000000000      0.5000000000"/>
  <vertex distance="3.714138460" upperboundary="0.8623187822" lowerboundary="-1.106211197" label="K" coord="0.3750000000      0.3750000000      0.7500000000"/>
  <vertex distance="3.714138460" upperboundary="0.8623187822" lowerboundary="-1.106211197" label="U" coord="0.6250000000      0.2500000000      0.6250000000"/>
  <vertex distance="3.930627631" upperboundary="0.8623187822" lowerboundary="-1.106211197" label="X" coord="0.5000000000       0.000000000      0.5000000000"/>
</bandstructure>
"""

