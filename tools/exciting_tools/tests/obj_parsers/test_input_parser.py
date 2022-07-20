"""
Test for the input.xml file parser
"""
import pytest

from excitingtools.exciting_obj_parsers.input_xml import parse_groundstate, parse_structure, parse_xs

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


def test_parse_groundstate_to_object():
    # just calls the function to see whether the initialization of object returns without errors
    ground_state = parse_groundstate(reference_input_str)
    assert ground_state.to_xml().tag == 'groundstate'


def test_parse_structure_to_object():
    # just calls the function to see whether the initialization of object returns without errors
    structure = parse_structure(reference_input_str)
    structure_xml = structure.to_xml()
    assert structure_xml.tag == 'structure'


def test_parse_xs_to_object():
    # just calls the function to see whether the initialization of object returns without errors
    xs = parse_xs(reference_input_str)
    xs_xml = xs.to_xml()
    assert xs_xml.tag == 'xs'


reference_input_str_without_xs = """<?xml version="1.0" encoding="UTF-8"?>
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

</input>
"""


def test_parse_xs_to_object_no_xs():
    xs = parse_xs(reference_input_str_without_xs)
    assert xs is None


reference_input_str_warning = """<?xml version="1.0" encoding="UTF-8"?>
<input>
  <xs xstype="BSE" 
   ngridq="3 3 3"
   vkloff="0.05 0.15 0.25" 
   nempty="1"
   broad="0.0073499"
   nosym="true">

   <transitions></transitions>
  </xs>

</input>
"""


def test_parse_xs_to_object_warning():
    with pytest.warns(UserWarning, match='Subelement transitions not yet supported. Its ignored...'):
        xs = parse_xs(reference_input_str_warning)
    xs_xml = xs.to_xml()
    assert xs_xml.tag == 'xs'
