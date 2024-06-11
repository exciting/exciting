"""Testing setup.diamond_phonon_g."""

import json
import os
from pathlib import Path

from excitingtools import ExcitingInputXML
from excitingtools.utils.test_utils import MockFile
from numpy.testing import assert_allclose

from excitingscripts.setup.diamond_phonon import point_specific_setup
from excitingscripts.setup.diamond_phonon_g import get_new_positions


def test_diamond_phonon_g(tmp_path: Path, input_xml_mock: MockFile) -> None:
    input_file = tmp_path / "input.xml"
    input_file.write_text(input_xml_mock.string)
    work_dir_name = "work_dir"

    os.chdir(tmp_path)
    point_specific_setup(get_new_positions)(0.5, 5, work_dir_name)

    work_dir = tmp_path / work_dir_name
    assert work_dir.exists()
    info_file = work_dir / "INFO-diamond-phonon.json"
    with open(info_file) as fid:
        info = json.load(fid)
    ref_info = {
        "Equilibrium lattice parameter (alat) in a.u.": 6.7274,
        "Maximum displacement (u,u,u); u in alat": 0.5,
        "Number of displacements": 5,
        "Volume of equilibrium unit cell in (a.u)^3": 76.11701721170597,
    }
    assert info == ref_info
    work_dir_content = {x.name for x in work_dir.iterdir()}
    ref_work_dir_content = {
        "INFO-diamond-phonon.json",
        "displ_0.25",
        "displ_0.5",
        "displ_1e-06",
        "displ_-0.25",
        "displ_-0.5",
    }
    assert work_dir_content == ref_work_dir_content

    input_obj = ExcitingInputXML.from_xml(work_dir / "displ_-0.5/input.xml")
    ref_positions = [[0.0, 0.0, 0.0], [-0.25, -0.25, -0.25]]
    assert_allclose(input_obj.structure.positions, ref_positions)

    input_obj = ExcitingInputXML.from_xml(work_dir / "displ_-0.25/input.xml")
    ref_positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    assert_allclose(input_obj.structure.positions, ref_positions)

    input_obj = ExcitingInputXML.from_xml(work_dir / "displ_1e-06/input.xml")
    ref_positions = [[0.0, 0.0, 0.0], [0.250001, 0.250001, 0.250001]]
    assert_allclose(input_obj.structure.positions, ref_positions)

    input_obj = ExcitingInputXML.from_xml(work_dir / "displ_0.25/input.xml")
    ref_positions = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    assert_allclose(input_obj.structure.positions, ref_positions)

    input_obj = ExcitingInputXML.from_xml(work_dir / "displ_0.5/input.xml")
    ref_positions = [[0.0, 0.0, 0.0], [0.75, 0.75, 0.75]]
    assert_allclose(input_obj.structure.positions, ref_positions)
