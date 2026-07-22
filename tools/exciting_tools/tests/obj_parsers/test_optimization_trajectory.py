"""Tests for the optimization trajectory parser."""

from pathlib import Path

from numpy.testing import assert_allclose

from excitingtools.exciting_obj_parsers.optimization_trajectory import parse_optimization_trajectory

ref_positions = [
    [
        [0.0, 0.0, 0.0],
        [0.31616, 0.31616, 0.62882],
        [0.74824, 0.74824, 0.10933],
        [0.56792, 0.56792, 0.51949],
        [0.83152, 0.83152, 0.87852],
        [0.48464, 0.48464, 0.7503],
        [0.82137, 0.82137, 0.42396],
        [0.49479, 0.49479, 0.20486],
        [0.15382, 0.15382, 0.57036],
        [0.16234, 0.16234, 0.05846],
    ],
    [
        [0.0, 0.0, 0.0],
        [0.31877553, 0.31877553, 0.636074],
        [0.74980014, 0.74980014, 0.11294626],
        [0.56897538, 0.56897538, 0.52312774],
        [0.83539467, 0.83539467, 0.87946348],
        [0.48338086, 0.48338086, 0.75661053],
        [0.82308023, 0.82308023, 0.42784111],
        [0.4956953, 0.4956953, 0.20823289],
        [0.15678175, 0.15678175, 0.57722036],
        [0.16199378, 0.16199378, 0.05885364],
    ],
]


def test_parse_optimization_trajectory(tmp_path: Path) -> None:
    info_out_str = (Path(__file__).parent / "INFO_opt.OUT").read_text()
    tmp_info_file = tmp_path / "INFO.OUT"
    tmp_info_file.write_text(info_out_str)

    trajectory = parse_optimization_trajectory(tmp_info_file)

    assert len(trajectory) == 2, "Expected 2 structures in the trajectory"
    assert_allclose(trajectory[0].positions, ref_positions[0])
    assert_allclose(trajectory[1].positions, ref_positions[1])
