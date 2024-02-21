import numpy as np
import pytest
from excitingscripts.execute.planar_average import execute_planar_average
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def vcl3d_xml_mock(tmp_path) -> MockFile:
    """ Mock 'VCL3D.xml' data.
    """
    vcl3d_xml_str = """
    <plot3d>
      <title>HCCF</title>
      <grid gridticks="5   5   4" origin="0.000       0.000       0.000" originrs="0.000       0.000       0.000">
        <axis name="a" label=" a" latexunit=" lattice coordinate" 
        graceunit=" graceunit" endpoint="1.000       0.000       0.000" delta="0.042       0.000       0.000"
        endpointrs="4.925       0.000       0.000"/>
        <axis name="b" label=" b" latexunit=" lattice coordinate" 
        graceunit=" graceunit" endpoint="0.000       1.000       0.000" delta="0.000       0.042       0.000" 
        endpointrs="2.462       4.265       0.000"/>
        <axis name="c" label=" c" latexunit=" lattice coordinate" 
        graceunit=" graceunit" endpoint="0.000       0.000       1.000" delta="0.000       0.000       0.008" 
        endpointrs="0.000       0.000      28.000"/>
        <value label=" Potential" latexunit=" E_h/(ea_0)" graceunit=" graceunit"/>
      </grid>
      <function n="78750">
        <row const="c" index="0">
          <row const="b" index="0"> -0.2010612988E-01   -0.2027721931E-01   -0.2053249465E-01   -0.2054464403E-01
             -0.2038953604E-01  </row>
          <row const="b" index="1"> -0.2027721931E-01   -0.2033858993E-01   -0.2032216053E-01   -0.2024525568E-01
             -0.2021821079E-01  </row>
          <row const="b" index="2"> -0.2053249465E-01   -0.2032216053E-01   -0.2009656614E-01   -0.2001183242E-01
             -0.2006196285E-01  </row>
          <row const="b" index="3"> -0.2054464403E-01   -0.2024525568E-01   -0.2001183242E-01   -0.1994198068E-01
             -0.1996753058E-01  </row>
          <row const="b" index="4"> -0.2038953604E-01   -0.2021821079E-01   -0.2006196285E-01   -0.1996753058E-01
             -0.1992338703E-01  </row>
        </row>
        <row const="c" index="1">
          <row const="b" index="0"> -0.2102700612E-01   -0.2107919668E-01   -0.2106905528E-01   -0.2082932286E-01
             -0.2053240504E-01  </row>
          <row const="b" index="1"> -0.2107919668E-01   -0.2097573818E-01   -0.2072680901E-01   -0.2045184283E-01
             -0.2033037757E-01  </row>
          <row const="b" index="2"> -0.2106905528E-01   -0.2072680901E-01   -0.2034535364E-01   -0.2013237449E-01
             -0.2012565658E-01  </row>
          <row const="b" index="3"> -0.2082932286E-01   -0.2045184283E-01   -0.2013237449E-01   -0.1998344517E-01
             -0.1996516084E-01  </row>
          <row const="b" index="4"> -0.2053240504E-01   -0.2033037757E-01   -0.2012565658E-01   -0.1996516084E-01
             -0.1987734695E-01  </row>
        </row>
        <row const="c" index="2">
          <row const="b" index="0"> -0.2126465084E-01   -0.2119527276E-01   -0.2094272461E-01   -0.2053443347E-01
             -0.2020066780E-01  </row>
          <row const="b" index="1"> -0.2119527276E-01   -0.2095793313E-01   -0.2054907830E-01   -0.2017164038E-01
             -0.2003464788E-01  </row>
          <row const="b" index="2"> -0.2094272461E-01   -0.2054907830E-01   -0.2011861398E-01   -0.1987385792E-01
             -0.1988518531E-01  </row>
          <row const="b" index="3"> -0.2053443347E-01   -0.2017164038E-01   -0.1987385792E-01   -0.1974621675E-01
             -0.1977379805E-01  </row>
          <row const="b" index="4"> -0.2020066780E-01   -0.2003464788E-01   -0.1988518531E-01   -0.1977379805E-01
             -0.1973230091E-01  </row>
        </row>
      </function>
    </plot3d>
    """

    vcl3d_xml_file = tmp_path / "VCL3D.xml"
    vcl3d_xml_file.write_text(vcl3d_xml_str)

    return MockFile(vcl3d_xml_file, vcl3d_xml_str)

def test_execute_planar_average(vcl3d_xml_mock, tmp_path):
    planar_average_data_ref = np.array([[0, -0.02027190],
                                        [14, -0.02068180],
                                        [28, -0.02027190]])

    assert np.allclose(execute_planar_average(f"{tmp_path}/VCL3D.xml", "z"), planar_average_data_ref)

