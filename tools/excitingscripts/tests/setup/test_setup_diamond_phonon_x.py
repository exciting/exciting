"""Testing setup.diamond_phonon_x."""

import json
import os
from pathlib import Path

from excitingtools import ExcitingInputXML
from excitingtools.utils.test_utils import MockFile
from numpy.testing import assert_allclose

from excitingscripts.setup.diamond_phonon import point_specific_setup
from excitingscripts.setup.diamond_phonon_x import get_new_positions, set_supercell

def test_diamond_phonon_x(tmp_path: Path, input_xml_mock: MockFile) -> None:
    input_file = tmp_path / "input.xml"
    input_file.write_text(input_xml_mock.string)
    work_dir_name = "work_dir"

    os.chdir(tmp_path)
    point_specific_setup(get_new_positions, set_supercell)(0.5, 5, work_dir_name, "TA")

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
        "X-phonon-calculation mode": "TA"
    }
    assert info == ref_info
    work_dir_content = {x.name for x in work_dir.iterdir()}
    ref_work_dir_content = {
        "INFO-diamond-phonon.json",
        "displ_0.25",
        "displ_0.5",
        "displ_1e-06",
    }
    assert work_dir_content == ref_work_dir_content

    input_obj = ExcitingInputXML.from_xml(work_dir / "displ_1e-06/input.xml")
    ref_positions = [[0.0, 1.4142135623730952e-06, 0.0], [0.0, 0.5000014142135624, 0.25],
                     [0.5, 0.49999858578643763, 0.5], [0.5, 0.9999985857864376, 0.75]]
    assert_allclose(input_obj.structure.positions, ref_positions)

    input_obj = ExcitingInputXML.from_xml(work_dir / "displ_0.25/input.xml")
    ref_positions = [[0.0, 0.3535533905932738, 0.0], [0.0, 0.8535533905932737, 0.25],
                     [0.5, 0.1464466094067262, 0.5], [0.5, 0.6464466094067263, 0.75]]
    assert_allclose(input_obj.structure.positions, ref_positions)

    input_obj = ExcitingInputXML.from_xml(work_dir / "displ_0.5/input.xml")
    ref_positions = [[0.0, 0.7071067811865476, 0.0], [0.0, 1.2071067811865475, 0.25],
                     [0.5, -0.20710678118654757, 0.5], [0.5,0.2928932188134524, 0.75]]
    assert_allclose(input_obj.structure.positions, ref_positions)
