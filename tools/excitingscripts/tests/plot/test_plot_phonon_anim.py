from pathlib import Path
from typing import Tuple

import pytest
from excitingscripts.plot.phonon_anim import parse_species_data
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def species_mock(tmp_path: Path) -> Tuple[MockFile, MockFile]:
    """Mock species.xml files."""
    c_xml_str = """<?xml version="1.0" encoding="UTF-8"?>
    <spdb xsi:noNamespaceSchemaLocation="../../xml/species.xsd" xmlns:xsi="...">
      <sp chemicalSymbol="C" name="carbon" z="-6.00000" mass="21894.16673">
        <muffinTin rmin="0.100000E-04" radius="1.4500" rinf="21.0932" radialmeshPoints="250"/>
        <atomicState n="1" l="0" kappa="1" occ="2.00000" core="true"/>
        <atomicState n="2" l="0" kappa="1" occ="2.00000" core="false"/>
        <atomicState n="2" l="1" kappa="1" occ="1.00000" core="false"/>
        <atomicState n="2" l="1" kappa="2" occ="1.00000" core="false"/>
        <basis>
          <default type="lapw" trialEnergy="0.1500" searchE="false"/>
          <custom l="0" type="apw+lo" trialEnergy="0.1500" searchE="true"/>
          <custom l="1" type="apw+lo" trialEnergy="0.1500" searchE="true"/>
        </basis>
      </sp>
    </spdb>"""

    o_xml_str = """<?xml version="1.0" encoding="UTF-8"?>
    <spdb xsi:noNamespaceSchemaLocation="../../xml/species.xsd" xmlns:xsi="...">
      <sp chemicalSymbol="O" name="oxygen" z="-8.00000" mass="29165.12203">
        <muffinTin rmin="0.100000E-04" radius="1.4500" rinf="17.0873" radialmeshPoints="250"/>
        <atomicState n="1" l="0" kappa="1" occ="2.00000" core="true"/>
        <atomicState n="2" l="0" kappa="1" occ="2.00000" core="false"/>
        <atomicState n="2" l="1" kappa="1" occ="2.00000" core="false"/>
        <atomicState n="2" l="1" kappa="2" occ="2.00000" core="false"/>
        <basis>
          <default type="lapw" trialEnergy="0.1500" searchE="false"/>
          <custom l="0" type="apw+lo" trialEnergy="0.1500" searchE="true"/>
          <custom l="1" type="apw+lo" trialEnergy="0.1500" searchE="true"/>
        </basis>
      </sp>
    </spdb>"""

    species_dir = tmp_path / "species"
    species_dir.mkdir()

    c_xml_file = species_dir / "C.xml"
    c_xml_file.write_text(c_xml_str)

    o_xml_file = species_dir / "O.xml"
    o_xml_file.write_text(o_xml_str)

    return MockFile(c_xml_file, c_xml_str), MockFile(o_xml_file, o_xml_str)


def test_parse_species_data(species_mock: Tuple[MockFile, MockFile], tmp_path: Path) -> None:
    unique_species = ["C", "O"]
    all_species = ["C", "C", "O"]

    species_data, natmax = parse_species_data(unique_species, tmp_path, all_species)

    assert len(species_data) == 2

    c_data = species_data[0]
    assert c_data["mass"] == 21894.16673
    assert c_data["atomic_number"] == 6
    assert c_data["count"] == 2

    o_data = species_data[1]
    assert o_data["mass"] == 29165.12203
    assert o_data["atomic_number"] == 8
    assert o_data["count"] == 1

    assert natmax == 2
