import pytest
from excitingscripts.setup.volume_optimization import setup_volume_optimization
from excitingtools.exciting_obj_parsers.input_xml import parse_input_xml
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def input_xml_mock(tmp_path) -> MockFile:
    """ Mock 'input.xml' data.
    """
    input_xml_str = """<?xml version="1.0" ?>
    <input>
        <title/>
        <structure speciespath="/home/exciting/species">
            <crystal scale="7.729">
                <basevect>0.5 0.5 0.0</basevect>
                <basevect>0.5 0.0 0.5</basevect>
                <basevect>0.0 0.5 0.5</basevect>
            </crystal>
            <species speciesfile="Ag.xml">
                <atom coord="0.0 0.0 0.0"> </atom>
            </species>
        </structure>
        <groundstate
            ngridk="8 8 8"
            rgkmax="7.5"
            swidth="0.01"
            xctype="GGA_PBE_SOL">
        </groundstate>

    </input>
    """

    input_xml_file = tmp_path / "input.xml"
    input_xml_file.write_text(input_xml_str)

    return MockFile(input_xml_file, input_xml_str)


def test_setup_volume_optimization(input_xml_mock, tmp_path):
    parsed_input = parse_input_xml(input_xml_mock.string)
    parsed_input.structure.crystal_properties.scale = 8.115450000000001

    setup_volume_optimization(input_xml_mock.full_path, 11, root_directory=tmp_path)
    parsed_input_volume_11 = parse_input_xml(tmp_path / "volume-11" / "input.xml")

    assert parsed_input.structure.crystal_properties.scale == parsed_input_volume_11.structure.crystal_properties.scale

