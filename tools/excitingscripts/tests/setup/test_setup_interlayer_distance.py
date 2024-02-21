import pytest
from excitingscripts.setup.interlayer_distance import setup_interlayer_distance
from excitingtools.exciting_obj_parsers.input_xml import parse_input_xml
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def input_xml_mock(tmp_path) -> MockFile:
    """ Mock 'input.xml' data.
    """
    input_xml_str = """<?xml version="1.0" ?>
    <input>
 
       <title>Graphite</title>
     
       <structure speciespath="/home/exciting/species">
          <crystal scale="1.88973" stretch="2.46 2.46 6.60">
             <basevect>0.5 0.8660254038 0.0</basevect>
             <basevect>1.0 0.0000000000 0.0</basevect>
             <basevect>0.0 0.0000000000 1.0</basevect>
          </crystal>
          <species speciesfile="C.xml" rmt="1.2">
             <atom coord="0.00000000 0.00000000 0.0"></atom>
             <atom coord="0.00000000 0.00000000 0.5"></atom>
             <atom coord="0.66666667 0.66666667 0.0"></atom>
             <atom coord="0.33333333 0.33333333 0.5"></atom>
          </species>
       </structure>
     
       <groundstate
           rgkmax="6.0"
           gmaxvr="20"
           ngridk="10 10 4"
           xctype="GGA_PBE"
           vdWcorrection="DFTD2">
       </groundstate>
     
    </input>
    """

    input_xml_file = tmp_path / "input.xml"
    input_xml_file.write_text(input_xml_str)

    return MockFile(input_xml_file, input_xml_str)


def test_setup_interlayer_distance(input_xml_mock, tmp_path):
    number_calculations = 4

    parsed_input = parse_input_xml(input_xml_mock.string)
    parsed_input.structure.lattice[2][2] = 1.2828512137937294

    setup_interlayer_distance(input_xml_mock.full_path, 5, 8, 4, 20,  root_directory=tmp_path)

    strain_values = []
    for i in range(number_calculations):
        with open(f"{tmp_path}/rundir-{i + 1}/strain-{i + 1}", "r") as f:
            strain_values.append(float(f.read()))

    parsed_input_rundir_4 = parse_input_xml(f"{tmp_path}/rundir-4/input.xml")

    assert strain_values == [5, 6, 7, 8]

    assert parsed_input.structure.lattice[2][2] == parsed_input_rundir_4.structure.lattice[2][2]

    with open(f"{tmp_path}/rundir-oo/strain-oo", "r") as f:
        assert float(f.read()) == 20

