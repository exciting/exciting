import numpy as np

from excitingtools.exciting_dict_parsers.properties_parser import parse_charge_density


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

def test_parse_charge_density():
    rho1 = parse_charge_density(RHO1_xml)
    ref = np.array([[0.00000000e+00, 1.98518837e+03],
                    [4.48811667e-02, 5.08882988e+02],
                    [8.97623334e-02, 1.51827164e+02],
                    [1.34643500e-01, 5.22636138e+01],
                    [1.79524667e-01, 2.43321622e+01],
                    [2.24405834e-01, 1.59499145e+01]])
    assert np.allclose(rho1, ref, atol=1.e-8)
