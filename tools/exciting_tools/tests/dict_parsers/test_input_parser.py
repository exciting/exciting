"""
Test for the input.xml file parser
"""
import pytest

from excitingtools.exciting_dict_parsers.input_parser import parse_element_xml, parse_structure, parse_input_xml

reference_input_str = """<?xml version="1.0" encoding="UTF-8"?>
<input sharedfs="true" 
  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" 
  xsi:noNamespaceSchemaLocation="https://xml.exciting-code.org/excitinginput.xsd">
  
  <title>Lithium Fluoride BSE</title>
  
  <structure speciespath="." autormt="false" epslat="1.0d-6">
    <crystal scale="1.0" stretch="1.0 1.0 1.0">
      <basevect>3.80402 3.80402 0.00000</basevect>
      <basevect>3.80402 0.00000 3.80402</basevect>
      <basevect>0.00000 3.80402 3.80402</basevect>
    </crystal>
    <species speciesfile="Li.xml" rmt="1.5">
      <atom coord="0.0000  0.0000  0.0000" bfcmt="0.0 0.0 0.0"/>
      <dfthalfparam 
        cut="3.90" 
        ampl="1" 
        exponent="8">
        <shell number="0" ionization="0.25" />
      </dfthalfparam>
    </species>
    <species speciesfile="F.xml">
      <atom coord="0.5000  0.5000  0.5000" lockxyz="false true false"/>
      <LDAplusU J="2.3" U="0.5" l="3"/>
    </species>
  </structure>
  
  <groundstate xctype="GGA_PBE" ngridk="4 4 4" epsengy="1d-7" outputlevel="high">
    <spin bfieldc="0 0 0" fixspin="total FSM"/>
    <OEP maxitoep="100"> </OEP>
  </groundstate>
  
  <properties>
    <dos 
      nsmdos="2"
      ngrdos="300"
      nwdos="1000"
      winddos="-0.3 0.3">
    </dos>
    <bandstructure>
      <plot1d>
        <path steps="100">
          <point coord="1.0     0.0     0.0" label="Gamma"/>
          <point coord="0.625   0.375   0.0" label="K"/>
          <point coord="0.5     0.5     0.0" label="X" breakafter="true"/>
          <point coord="0.0     0.0     0.0" label="Gamma"/>
          <point coord="0.5     0.0     0.0" label="L"/>
        </path>
      </plot1d>
    </bandstructure>
  </properties>

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


def test_parse_title():
    assert parse_element_xml(reference_input_str, tag="title") == "Lithium Fluoride BSE"


def test_parse_groundstate():
    ground_state = parse_element_xml(reference_input_str, tag="groundstate")
    assert ground_state == {
        'xctype': 'GGA_PBE', 'ngridk': [4, 4, 4],
        'epsengy': 1e-7, 'outputlevel': 'high',
        'spin': {'bfieldc': [0, 0, 0], 'fixspin': 'total FSM'},
        'OEP': {'maxitoep': 100}
    }


def test_parse_groundstate_from_gs_root():
    ground_state = parse_element_xml('<groundstate xctype="GGA_PBE" ngridk="4 4 4" epsengy="1d-7" outputlevel="high"/>',
                                     tag="groundstate")
    assert ground_state == {
        'xctype': 'GGA_PBE', 'ngridk': [4, 4, 4],
        'epsengy': 1e-7, 'outputlevel': 'high'
    }


def test_parse_structure():
    structure = parse_structure(reference_input_str)
    structure_ref = {
        'atoms': [{'species': 'Li', 'position': [0.0, 0.0, 0.0],
                   'bfcmt': [0.0, 0.0, 0.0]},
                  {'species': 'F', 'position': [0.5, 0.5, 0.5],
                   'lockxyz': [False, True, False]}],
        'lattice': [[3.80402, 3.80402, 0.0],
                    [3.80402, 0.0, 3.80402],
                    [0.0, 3.80402, 3.80402]],
        'species_path': '.',
        'crystal_properties': {'scale': 1.0, 'stretch': [1.0, 1.0, 1.0]},
        'species_properties': {'F': {'LDAplusU': {'J': 2.3, 'U': 0.5, 'l': 3}},
                               'Li': {'dfthalfparam': {'ampl': 1,
                                                       'cut': 3.9,
                                                       'exponent': 8,
                                                       'shell': [{'ionization': 0.25,
                                                                 'number': 0}]},
                                      'rmt': 1.5}},
        'autormt': False,
        'epslat': 1.0e-6,
    }
    assert structure_ref == structure


def test_parse_xs():
    xs = parse_element_xml(reference_input_str, tag="xs")
    xs_ref = {
        'xstype': 'BSE',
        'ngridq': [3, 3, 3],
        'vkloff': [0.05, 0.15, 0.25],
        'nempty': 1,
        'broad': 0.0073499,
        'nosym': True,
        'energywindow': {'intv': [0.0, 1.0], 'points': 50},
        'screening': {'screentype': 'full', 'nempty': 115},
        'BSE': {'bsetype': 'singlet', 'nstlbse': [1, 5, 1, 2], 'aresbse': False},
        'qpointset': [[0.0, 0.0, 0.0]],
        'plan': ['screen', 'bse']
    }
    assert xs_ref == xs
    assert isinstance(xs["ngridq"][0], int)


input_ref_parsed_keys = {'title', 'groundstate', 'structure', 'xs', 'sharedfs', "properties"}


def test_parse_input_xml():
    parsed_data = parse_element_xml(reference_input_str)
    assert set(parsed_data.keys()) == input_ref_parsed_keys
    assert parsed_data['sharedfs']


def test_parse_input_xml_directly():
    parsed_data = parse_input_xml(reference_input_str)
    assert set(parsed_data.keys()) == input_ref_parsed_keys


def test_parse_missing_tag():
    with pytest.raises(ValueError, match="Your specified input has no tag missing_tag"):
        parse_element_xml(reference_input_str, "missing_tag")


def test_parse_input_xml_with_tag():
    parsed_data = parse_element_xml(reference_input_str, tag="input")
    assert set(parsed_data.keys()) == input_ref_parsed_keys


def test_parse_properties():
    properties = parse_element_xml(reference_input_str, tag="properties")
    properties_ref = {
        'dos': {"nsmdos": 2, "ngrdos": 300, "nwdos": 1000, "winddos": [-0.3, 0.3]},
        'bandstructure': {"plot1d": {"path": {
            "steps": 100,
            "point":
                [
                    {"coord": [1, 0, 0], "label": "Gamma"},
                    {"coord": [0.625, 0.375, 0], "label": "K"},
                    {"coord": [0.5, 0.5, 0], "label": "X", "breakafter": True},
                    {"coord": [0, 0, 0], "label": "Gamma"},
                    {"coord": [0.5, 0, 0], "label": "L"}
                ]}}}
    }
    assert properties_ref == properties
