import pytest
from excitingscripts.setup.dos_band_structure import setup_dos_band_structure

from excitingtools.exciting_dict_parsers.input_parser import parse_element_xml
from excitingtools.exciting_obj_parsers.input_xml import parse_input_xml
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def input_xml_mock(tmp_path) -> MockFile:
    """ Mock 'input.xml' data.
    """
    input_xml_str = """<?xml version="1.0" ?>
    <input>

       <title>Diamond</title>
       
	   <structure speciespath="/home/exciting/species">
		  <crystal scale="6.7274">
			 <basevect>0.0 0.5 0.5</basevect>
			 <basevect>0.5 0.0 0.5</basevect>
			 <basevect>0.5 0.5 0.0</basevect>
		  </crystal>
		  <species speciesfile="C.xml">
			 <atom coord="0.0 0.0 0.0"> </atom>
			 <atom coord="0.25 0.25 0.25"> </atom>
		  </species>
	   </structure>
	   
       <groundstate
           ngridk="4 4 4"
           outputlevel="normal"
           xctype="GGA_PBE_SOL">
       </groundstate>

    </input>
    """

    input_xml_file = tmp_path / "input.xml"
    input_xml_file.write_text(input_xml_str)

    return MockFile(input_xml_file, input_xml_str)


def test_setup_dos_band_structure(input_xml_mock, tmp_path):
    properties_ref = {
        "dos": {"nsmdos": 2, "nwdos": 1000, "winddos": [-0.3, 0.3]},
        "bandstructure": {"plot1d": {"path": {
            "steps": 100,
            "point":
                [
                    {"coord": [0.0, 0.0, 0.0], "label": "G"},
                    {"coord": [0.5, 0.0, 0.5], "label": "X"},
                    {"coord": [0.5, 0.25, 0.75], "label": "W"},
                    {"coord": [0.375, 0.375, 0.75], "label": "K"},
                    {"coord": [0.0, 0.0, 0.0], "label": "G"},
                    {"coord": [0.5, 0.5, 0.5], "label": "L"},
                    {"coord": [0.625, 0.25, 0.625], "label": "U"},
                    {"coord": [0.5, 0.25, 0.75], "label": "W"},
                    {"coord": [0.5, 0.5, 0.5], "label": "L"},
                    {"coord": [0.375, 0.375, 0.75], "label": "K", "breakafter": True},
                    {"coord": [0.625, 0.25, 0.625], "label": "U"},
                    {"coord": [0.5, 0.0, 0.5], "label": "X"}]}}}}

    setup_dos_band_structure(input_xml_mock.full_path, root_directory=tmp_path)
    parsed_input = parse_input_xml(tmp_path / "input.xml")
    parsed_input_properties = parsed_input.properties.to_xml()  # pylint: disable=no-member

    assert properties_ref == parse_element_xml(parsed_input_properties)
