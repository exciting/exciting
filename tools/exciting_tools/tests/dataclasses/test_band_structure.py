# TODO(Mara) Remove call to parse_band_structure and mock inputs to BandData using np, list, etc.
from excitingtools.exciting_obj_parsers.ks_band_structure import parse_band_structure
from excitingtools.dataclasses.band_structure import BandData

import numpy as np


def test_class_band_data_xticks_and_labels():

    xticks_ref = [0.000000000, 0.6123238446, 0.9184857670, 1.134974938, 1.784442453, 2.314730457,
                    2.689700702, 2.906189873, 3.339168216, 3.714138460, 3.930627631]
    unicode_gamma = '\u0393'
    labels_ref = [unicode_gamma, "X", "W", "K", unicode_gamma, "L", "U", "W", "L", "K,U", "X"]

    band_data: BandData = parse_band_structure(band_structure_xml)
    xticks, labels = band_data.band_path()

    assert np.allclose(xticks, xticks_ref)
    assert labels == labels_ref


# TODO(Mara) Remove after mocking data from it
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
