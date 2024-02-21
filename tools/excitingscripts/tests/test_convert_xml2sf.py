import filecmp

import pytest
from excitingscripts.convert_xml2xsf import convert_xml2xsf
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def WF3D_xml_mock(tmp_path) -> MockFile:
    """ Mock 'WF3D.xml' data.
    """
    WF3D_xml_str = """<?xml version="1.0" ?>
    <plot3d>
      <title>Bulk Silicon: Plot example</title>
      <grid gridticks="3    3    3" origin="-0.500      -0.500      -0.500" originrs="-5.130      -5.130      -5.130">
        <axis name="a" label=" a" latexunit=" lattice coordinate" graceunit=" graceunit" endpoint="1.000       0.000       0.000" delta="0.050       0.000       0.000" endpointrs="0.000       5.130       5.130"/>
        <axis name="b" label=" b" latexunit=" lattice coordinate" graceunit=" graceunit" endpoint="0.000       1.000       0.000" delta="0.000       0.050       0.000" endpointrs="5.130       0.000       5.130"/>
        <axis name="c" label=" c" latexunit=" lattice coordinate" graceunit=" graceunit" endpoint="0.000       0.000       1.000" delta="0.000       0.000       0.050" endpointrs="5.130       5.130       0.000"/>
        <value label=" Wave Function Norm Squared" latexunit=" " graceunit=" graceunit"/>
      </grid>
      <function n="9261">
        <row const="c" index="0">
          <row const="b" index="0">  0.1694892838E-19    0.2627258983E-07    0.1372458919E-06    0.3875770767E-06  </row>
          <row const="b" index="1">  0.2627262588E-07    0.1522766640E-06    0.4176435376E-06    0.8549474969E-06  </row>
          <row const="b" index="2">  0.1372458885E-06    0.4176433737E-06    0.8586449835E-06    0.1434067503E-05  </row>
          <row const="b" index="3">  0.3875765866E-06    0.8549465150E-06    0.1434066537E-05    0.2013978713E-05  </row>
        </row>
    
      </function>
    </plot3d>
    """

    input_xml_file = tmp_path / "WF3D.xml"
    input_xml_file.write_text(WF3D_xml_str)

    return MockFile(input_xml_file, WF3D_xml_str)

def test_convert_xml2xsf(WF3D_xml_mock, tmp_path):
    WF3D_xsf_ref_str = """
BEGIN_BLOCK_DATAGRID_3D
Bulk_Silicon:_Plot_example      
BEGIN_DATAGRID_3D
3    3    3
-2.71467801 -2.71467801 -2.71467801 
0 2.71467801 2.71467801 
2.71467801 0 2.71467801 
2.71467801 2.71467801 0 
  0.1694892838E-19    0.2627258983E-07    0.1372458919E-06    0.3875770767E-06  
  0.2627262588E-07    0.1522766640E-06    0.4176435376E-06    0.8549474969E-06  
  0.1372458885E-06    0.4176433737E-06    0.8586449835E-06    0.1434067503E-05  
  0.3875765866E-06    0.8549465150E-06    0.1434066537E-05    0.2013978713E-05  
END_DATAGRID_3D
END_BLOCK_DATAGRID_3D
"""

    with open(tmp_path / "WF3D_ref.xsf", "w") as fid:
        fid.write(WF3D_xsf_ref_str)

    convert_xml2xsf(WF3D_xml_mock.full_path, dimension="3D")

    assert filecmp.cmp(tmp_path / "WF3D_ref.xsf", tmp_path / "WF3D.xsf")
