"""Test parser factory."""

from pathlib import Path

from numpy.testing import assert_allclose

from excitingtools.exciting_dict_parsers.parser_factory import parse


def test_parse(tmp_path: Path) -> None:
    file = tmp_path / "EPSILON_11.OUT"
    file.write_text("a\n0 1 0\n2 0 0")
    parsed_data = parse(file.as_posix())
    assert set(parsed_data) == {"energy", "im", "re"}
    assert_allclose(parsed_data["energy"], [0.0, 2.0])
    assert_allclose(parsed_data["im"], [0.0, 0.0])
    assert_allclose(parsed_data["re"], [1.0, 0.0])

    file = tmp_path / "EPSILON_BSE-NAR_TDA-BAR_OC11.OUT"
    file.write_text("a\n" * 14 + "0 0 0 1\n0 0 0 2")
    parsed_data = parse(file.as_posix())
    assert set(parsed_data) == {
        "frequency",
        "imag_oscillator_strength",
        "real_oscillator_strength",
        "real_oscillator_strength_kkt",
    }
    assert_allclose(parsed_data["frequency"], [0.0, 0.0])
    assert_allclose(parsed_data["imag_oscillator_strength"], [0.0, 0.0])
    assert_allclose(parsed_data["real_oscillator_strength"], [0.0, 0.0])
    assert_allclose(parsed_data["real_oscillator_strength_kkt"], [1.0, 2.0])

    file = tmp_path / "Z_11.OUT"
    file.write_text("0 0 1 2\n0 0 3 4")
    parsed_data = parse(file.as_posix())
    assert set(parsed_data) == {"im", "mu", "re", "temperature"}
    assert_allclose(parsed_data["im"], [2.0, 4.0])
    assert_allclose(parsed_data["mu"], [0.0, 0.0])
    assert_allclose(parsed_data["re"], [1.0, 3.0])
    assert_allclose(parsed_data["temperature"], [0.0, 0.0])
