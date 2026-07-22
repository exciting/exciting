import os

import numpy as np
import pytest
from excitingscripts.setup.convergence_test import setup_convergence_test
from excitingtools.exciting_obj_parsers.input_xml import parse_input_xml
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def input_xml_mock(tmp_path) -> MockFile:
    """ Mock 'input.xml' data.
    """
    input_xml_4_5_str = """<?xml version="1.0" encoding="UTF-8"?>
    <input>
        <title/>
        <structure speciespath="/home/exciting/species">
            <crystal scale="7.7201">
                <basevect>0.5 0.5 0.0</basevect>
                <basevect>0.5 0.0 0.5</basevect>
                <basevect>0.0 0.5 0.5</basevect>
            </crystal>
            <species speciesfile="Ag.xml">
                <atom coord="0.0 0.0 0.0"> </atom>
            </species>
        </structure>
        <groundstate
           xctype="GGA_PBE_SOL"
           ngridk="4 4 4"
           rgkmax="5">
        </groundstate>

    </input>
    """

    os.makedirs(os.path.dirname(tmp_path / "4_5/input.xml"), exist_ok=True)

    input_xml_4_5_file = tmp_path / "4_5/input.xml"
    input_xml_4_5_file.write_text(input_xml_4_5_str)

    return MockFile(input_xml_4_5_file, input_xml_4_5_str)


def test_setup_convergence_test(input_xml_mock, tmp_path):
    parsed_input = parse_input_xml(input_xml_mock.string)
    parsed_input.groundstate.ngridk = [6, 6, 6]
    parsed_input.groundstate.rgkmax = 6

    setup_convergence_test(input_xml_mock.full_path, 4, 6, 2, 5, 6, 1, root_directory=tmp_path)
    parsed_input_6_6 = parse_input_xml(f"{tmp_path}/6_6/input.xml")

    assert np.allclose(parsed_input.groundstate.ngridk, parsed_input_6_6.groundstate.ngridk)
    assert parsed_input.groundstate.rgkmax == parsed_input_6_6.groundstate.rgkmax
