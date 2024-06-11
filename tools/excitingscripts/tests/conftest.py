"""General re-used pytest fixtures."""
import pytest
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
            <basevect>0.5 0.5 0.0</basevect>
             <basevect>0.5 0.0 0.5</basevect>
             <basevect>0.0 0.5 0.5</basevect>
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

@pytest.fixture
def phonon_results_mock(tmp_path) -> MockFile:
    """ Mock "phonon_results.json" data.
    """
    phonon_results_str = """
    {
        "exciting_version": "/home/exciting/bin/exciting_smp",
        "results": {
            "-0.025": {
                "energy": -75.9459459515,
                "force": 0.1244712479
            },
            "-0.02": {
                "energy": -75.9570243257,
                "force": 0.0943004005
            },
            "-0.015": {
                "energy": -75.9651916303,
                "force": 0.0670517691
            },
            "-0.01": {
                "energy": -75.9707282293,
                "force": 0.0424301339
            },
            "-0.005": {
                "energy": -75.973885787,
                "force": 0.020170169
            },
            "1e-06": {
                "energy": -75.9748914348,
                "force": 2.78895e-05
            },
            "0.005": {
                "energy": -75.9739487087,
                "force": -0.0182042428
            },
            "0.01": {
                "energy": -75.9712392101,
                "force": -0.0347361639
            },
            "0.015": {
                "energy": -75.9669252271,
                "force": -0.049740614
            },
            "0.02": {
                "energy": -75.9611524321,
                "force": -0.0633704292
            },
            "0.025": {
                "energy": -75.9540536417,
                "force": -0.0757527927
            }
        }
    }
    """

    phonon_results_file = tmp_path / "phonon_results.json"
    phonon_results_file.write_text(phonon_results_str)

    return MockFile(phonon_results_file, phonon_results_str)

@pytest.fixture
def info_diamond_phonon_mock(tmp_path) -> MockFile:
    """ Mock "INFO-diamond-phonon.json" data.
    """
    info_diamond_phonon_str = """
    {
        "Equilibrium lattice parameter (alat) in a.u.": 6.7468,
        "Maximum displacement (u,u,u); u in alat": 0.025,
        "Number of displacements": 11,
        "Volume of equilibrium unit cell in (a.u)^3": 76.77742058180796
    }
    """

    info_diamond_phonon_file = tmp_path / "INFO-diamond-phonon.json"
    info_diamond_phonon_file.write_text(info_diamond_phonon_str)

    return MockFile(info_diamond_phonon_file, info_diamond_phonon_str)
