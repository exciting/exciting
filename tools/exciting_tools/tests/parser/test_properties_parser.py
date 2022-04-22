import numpy as np
import pytest

from excitingtools.parser.properties_parser import parse_band_structure_to_arrays, parse_charge_density


def test_parse_band_structure_to_arrays():
    band_data = parse_band_structure_to_arrays(band_structure_xml)

    assert band_data.k_points.shape[0] == band_data.bands.shape[0], (
        "First dim of bands array equals the number of k-sampling points in the band structure")
    assert band_data.k_points.shape[0] == 6, "sampling points per band"
    assert band_data.n_bands == 2, "band_structure_xml contains two bands"

    ref_k_points = [0.,  0.04082159, 0.08164318, 0.12246477, 0.16328636, 0.20410795]

    ref_bands = np.array([[-0.45003454, -0.00937631],
                          [-0.44931675, -0.01419609],
                          [-0.44716535, -0.02681183],
                          [-0.44358550, -0.04401707],
                          [-0.43858583, -0.06361773],
                          [-0.43217908, -0.0844593]])

    assert np.allclose(band_data.k_points, ref_k_points, atol=1.e-8)
    assert np.allclose(band_data.bands, ref_bands, atol=1.e-8)


def test_parse_band_structure_to_arrays_vertices():
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

    band_data = parse_band_structure_to_arrays(band_structure_xml)
    assert band_data.vertices == vertices_ref

    with pytest.raises(NotImplementedError) as error:
        band_data.band_path()

    assert error.value.args[0] == (
        "Getting a list of high symmetry points from self.vertices not implemented")


def test_parse_charge_density():
    rho1 = parse_charge_density(RHO1_xml)
    ref = np.array([[0.00000000e+00, 1.98518837e+03],
                    [4.48811667e-02, 5.08882988e+02],
                    [8.97623334e-02, 1.51827164e+02],
                    [1.34643500e-01, 5.22636138e+01],
                    [1.79524667e-01, 2.43321622e+01],
                    [2.24405834e-01, 1.59499145e+01]])
    assert np.allclose(rho1, ref, atol=1.e-8)


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

# Most of the values have been removed to shorten the test

RHO1_xml = """<?xml version="1.0" encoding="UTF-8"?>
<plot1d>
  <title>silicon-primitive-PBEsol</title>
  <grid>
    <axis label=" Distance" latexunit=" a_0" graceunit=" graceunit"/>
    <axis label=" Density" latexunit=" ???" graceunit=" graceunit"/>
    <function name="">
      <point distance="0.000000000" value="1985.188367"/>
      <point distance="0.4488116670E-01" value="508.8829884"/>
      <point distance="0.8976233341E-01" value="151.8271641"/>
      <point distance="0.1346435001" value="52.26361379"/>
      <point distance="0.1795246668" value="24.33216220"/>
      <point distance="0.2244058335" value="15.94991449"/>
    </function>
    <vertex distance="0.000000000" upperboundary="1985.188367" lowerboundary="0.8456287531E-01" label="" coord="0.000000000       0.000000000       0.000000000"/>
    <vertex distance="4.443235504" upperboundary="1985.188367" lowerboundary="0.8456287531E-01" label="" coord="0.2500000000      0.2500000000      0.2500000000"/>
  </grid>
</plot1d>
"""
