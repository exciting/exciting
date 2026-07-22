import pytest
from excitingscripts.setup.dft_05 import setup_dft_05
from excitingtools.exciting_obj_parsers.input_xml import parse_input_xml
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def input_xml_mock(tmp_path) -> MockFile:
    """ Mock 'input.xml' data.
    """
    input_xml_str = """<?xml version="1.0" ?>
    <input>

        <title>Bulk Silicon: LDA-1/2</title>
       
        <structure speciespath="/home/exciting/species">
           <crystal scale="10.26">
             <basevect>0.0   0.5   0.5</basevect>
             <basevect>0.5   0.0   0.5</basevect>
             <basevect>0.5   0.5   0.0</basevect>
          </crystal>
          <species speciesfile="Si.xml" rmt="2.1">
             <atom coord="0.00  0.00  0.00"></atom>
             <atom coord="0.25  0.25  0.25"></atom>
             <dfthalfparam>
                <shell cut="3.90" ampl="1" exponent="8" number="0" ionization="0.25" />
             </dfthalfparam>
          </species>
       </structure>
     
       <groundstate
          do="fromscratch"
          rgkmax="7.0"
          gmaxvr="14"
          ngridk="6 6 6"
          outputlevel="high"
          xctype="LDA_PW">
          <dfthalf printVSfile="false"/>
       </groundstate>
     
    </input>
    """

    input_xml_file = tmp_path / "input.xml"
    input_xml_file.write_text(input_xml_str)

    return MockFile(input_xml_file, input_xml_str)


def test_setup_dft_05(input_xml_mock, tmp_path):
    number_calculations = 4

    parsed_input = parse_input_xml(input_xml_mock.string)
    species = parsed_input.structure.species[0]

    setup_dft_05(input_xml_mock.full_path, 0, 3, 4, species, root_directory=tmp_path)

    strain_values = []
    for i in range(number_calculations):
        with open(f"{tmp_path}/rundir-{i + 1}/strain-{i + 1}", "r") as f:
            strain_values.append(float(f.read()))

    parsed_input_rundir_4 = parse_input_xml(f"{tmp_path}/rundir-4/input.xml")

    assert strain_values == [0, 1, 2, 3]

    shell = parsed_input_rundir_4.structure.species_properties[species].dfthalfparam.shell[0]
    assert shell.cut == strain_values[3]
