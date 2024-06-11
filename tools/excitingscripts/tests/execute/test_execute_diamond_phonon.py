import os
from typing import Tuple

import excitingscripts
import pytest
from excitingscripts.execute.diamond_phonon import execute_diamond_phonon
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def info_xml_mock(tmp_path) -> Tuple[MockFile, MockFile, MockFile, MockFile, MockFile]:
    """ Mock 'info.xml' data.
    """
    info_xml_1_str = """
    <info>
      <groundstate status="finished">
        <scl>
          <iter iteration="1">
            <energies totalEnergy="-75.9459459371"></energies>
          </iter>
          <iter iteration="2">
            <energies totalEnergy="-75.9459459515"></energies>
          </iter>
          <structure>
            <species chemicalSymbol="C">
              <atom x="0.00000000000" y="0.00000000000" z="0.00000000000">
                <forces Magnitude="0.163297799372E-01">
                  <totalforce x="-0.124471247945" y="-0.124471247945" z="-0.124471247945"/>
                </forces>
              </atom>
              <atom x="0.00000000000" y="0.500000000000" z="0.250000000000">
                <forces Magnitude="0.163297799372E-01">
                  <totalforce x="0.124471247945" y="0.124471247945" z="0.124471247945"/>
                </forces>
              </atom>
            </species>
          </structure>
        </scl>
      </groundstate>
    </info>  
    """
    info_xml_2_str = info_xml_1_str.replace("-75.9459459515", "-75.9682724533")
    info_xml_2_str = info_xml_2_str.replace("0.124471247945", "0.544298359284E-01")

    info_xml_3_str = info_xml_1_str.replace("-75.9459459515", "-75.9540536417")
    info_xml_3_str = info_xml_3_str.replace("0.124471247945", "-0.757527927106E-01")

    info_xml_4_str = info_xml_1_str.replace("-75.9459459515", "-75.9692733329")
    info_xml_4_str = info_xml_4_str.replace("0.124471247945", "-0.424191910181E-01")

    info_xml_5_str = info_xml_1_str.replace("-75.9459459515", "-75.9748914348")
    info_xml_5_str = info_xml_5_str.replace("0.124471247945", "0.278894792514E-04")

    run_directories = ["displ_-0.025", "displ_-0.0125", "displ_0.025", "displ_0.0125", "displ_1e-06"]
    for directory in run_directories:
        os.makedirs(os.path.dirname(tmp_path / f"{directory}/info.xml"), exist_ok=True)


    info_xml_1_file = tmp_path / "displ_-0.025/info.xml"
    info_xml_1_file.write_text(info_xml_1_str)

    info_xml_2_file = tmp_path / "displ_-0.0125/info.xml"
    info_xml_2_file.write_text(info_xml_2_str)

    info_xml_3_file = tmp_path / "displ_0.025/info.xml"
    info_xml_3_file.write_text(info_xml_3_str)

    info_xml_4_file = tmp_path / "displ_0.0125/info.xml"
    info_xml_4_file.write_text(info_xml_4_str)

    info_xml_5_file = tmp_path / "displ_1e-06/info.xml"
    info_xml_5_file.write_text(info_xml_5_str)


    return MockFile(info_xml_1_file, info_xml_1_str), MockFile(info_xml_2_file, info_xml_2_str), \
        MockFile(info_xml_3_file, info_xml_3_str), MockFile(info_xml_4_file, info_xml_4_str), \
        MockFile(info_xml_5_file, info_xml_5_str)


def test_execute_diamond_phonon(monkeypatch, info_xml_mock, info_diamond_phonon_mock, tmp_path):

    def skip_run_exciting(*args, **kwargs):
        pass

    monkeypatch.setattr(excitingscripts.execute.diamond_phonon, 'run_exciting', skip_run_exciting)

    results_dict = execute_diamond_phonon(work_dir=tmp_path)["results"]

    results_dict_ref = {-0.025: {'energy': -75.9459459515, 'force': 0.1244712479},
                        -0.0125: {'energy': -75.9682724533, 'force': 0.0544298359},
                        1e-06: {'energy': -75.9748914348, 'force': 2.78895e-05},
                        0.0125: {'energy': -75.9692733329, 'force': -0.042419191},
                        0.025: {'energy': -75.9540536417, 'force': -0.0757527927}}

    assert results_dict == results_dict_ref
