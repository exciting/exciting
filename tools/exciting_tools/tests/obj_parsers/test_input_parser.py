"""
Test for the input.xml file parser
"""

from excitingtools.exciting_obj_parsers.input_xml import parse_input_xml

reference_input_str = """<?xml version="1.0" encoding="UTF-8"?>
<input>
  
  <title>Lithium Fluoride BSE</title>
  
  <structure speciespath="." autormt="false" epslat="1.0d-6">
    <crystal scale="1.0" stretch="1.0">
      <basevect>3.80402 3.80402 0.00000</basevect>
      <basevect>3.80402 0.00000 3.80402</basevect>
      <basevect>0.00000 3.80402 3.80402</basevect>
    </crystal>
    <species speciesfile="Li.xml" rmt="1.5">
      <atom coord="0.0000  0.0000  0.0000" bfcmt="0.0 0.0 0.0"/>
    </species>
    <species speciesfile="F.xml">
      <atom coord="0.5000  0.5000  0.5000" lockxyz="false true false"/>
    </species>
  </structure>
  
  <groundstate xctype="GGA_PBE" ngridk="4 4 4" epsengy="1d-7" outputlevel="high"/>

  <xs xstype="BSE" 
   ngridq="3 3 3"
   vkloff="0.05 0.15 0.25" 
   nempty="1"
   broad="0.0073499"
   nosym="true">

   <energywindow intv="0.0 1.0" 
    points="50" />

   <screening screentype="full"
    nempty="115" />

   <BSE bsetype="singlet"
    nstlbse="1 5 1 2" 
    aresbse="false"/>

   <qpointset>
      <qpoint>0.0 0.0 0.0</qpoint>
   </qpointset>
   
   <plan>
      <doonly task="screen" />
      <doonly task="bse" />
   </plan>
  </xs>
  
</input>
"""


def test_parse_input_xml_to_object():
    input_xml = parse_input_xml(reference_input_str)
    assert set(vars(input_xml)) == {'xs', 'groundstate', 'structure', 'title'}
    assert input_xml.to_xml_str().startswith('<?xml version="1.0" ?>\n<input>\n\t<title>Lithium Fluoride BSE')
