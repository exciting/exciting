"""
Test for the RT_TDDFT_parser
"""

import numpy as np
import pytest

from excitingtools.exciting_dict_parsers.RT_TDDFT_parser import (
    parse_atom_position_velocity_force,
    parse_eigval_screenshots,
    parse_etot,
    parse_force,
    parse_jind,
    parse_nexc,
    parse_occupations,
    parse_proj_screenshots,
    parse_rttddft_polarization,
)

atom_str = """0.00  0.00  0.00  0.00  0.00  0.00  0.00 -0.00 -0.02 -0.03
0.20 -0.01 -0.01 -0.00 -0.13 -0.03 -0.01 -0.65 -0.02 -0.42
0.40 -0.05 -0.04 -0.00 -0.22 -0.05 -0.03 -0.63 -0.02 -0.45
0.60 -0.09 -0.09 -0.01 -0.27 -0.02 -0.04 -0.06 -0.02 -0.45
0.80 -0.16 -0.16 -0.02 -0.33 -0.06 -0.06 -0.96 -0.03 -0.45
1.00 -0.23 -0.26 -0.04 -0.41 -0.05 -0.08 -0.53 -0.03 -0.44
"""


atom_ref = {
    "Time": np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]),
    "x": np.array([0.0, -0.01, -0.05, -0.09, -0.16, -0.23]),
    "y": np.array([0.0, -0.01, -0.04, -0.09, -0.16, -0.26]),
    "z": np.array([0.0, 0.0, 0.0, -0.01, -0.02, -0.04]),
    "vx": np.array([0.0, -0.13, -0.22, -0.27, -0.33, -0.41]),
    "vy": np.array([0.0, -0.03, -0.05, -0.02, -0.06, -0.05]),
    "vz": np.array([0.0, -0.01, -0.03, -0.04, -0.06, -0.08]),
    "Fx": np.array([0.0, -0.65, -0.63, -0.06, -0.96, -0.53]),
    "Fy": np.array([-0.02, -0.02, -0.02, -0.02, -0.03, -0.03]),
    "Fz": np.array([-0.03, -0.42, -0.45, -0.45, -0.45, -0.44]),
}


def test_parse_atom_position_velocity_force(tmp_path) -> None:
    file_path = tmp_path / "ATOM_0001.OUT"
    file_path.write_text(atom_str)
    atom = parse_atom_position_velocity_force(file_path.as_posix())
    keys = atom.keys()
    assert keys == atom_ref.keys()
    assert all([np.allclose(atom[key], atom_ref[key]) for key in keys])


eigval_str = """ik =       1
      1     -5.969546700602
      2     -0.241574097700
      3      0.020505678329
      4      0.174433598397
      5      0.182344505943
      6      0.550215820881
      7      0.720938702290
      8      1.094268190606
      9      1.112731400960
     10      1.187472826613

ik =       2
      1     -5.967597065501
      2     -0.280574330557
      3     -0.005696262302
      4      0.319065586764
      5      0.330768479700
      6      0.741425584593

"""


eigs_ref = {
    "kpoints": [
        {
            "ik": 1,
            "eigenvalues": [
                -5.969546700602,
                -0.241574097700,
                0.020505678329,
                0.174433598397,
                0.182344505943,
                0.550215820881,
                0.720938702290,
                1.094268190606,
                1.112731400960,
                1.187472826613,
            ],
        },
        {
            "ik": 2,
            "eigenvalues": [
                -5.967597065501,
                -0.280574330557,
                -0.005696262302,
                0.319065586764,
                0.330768479700,
                0.741425584593,
            ],
        },
    ]
}


def test_parse_eigval_screenshots(tmp_path) -> None:
    file_path = tmp_path / "EIGVAL_20.OUT"
    file_path.write_text(eigval_str)
    eigs = parse_eigval_screenshots(file_path.as_posix())
    assert eigs == eigs_ref


etot_str = """ Time  ETOT Madelung Eigenvalues-Core Eigenvalues-Valence Exchange Correlation XC-potential Coulomb pot. energy
0.0 -9.7 -9.8 -2.3 1.2 -0.3 -0.2 1.3 0.4
0.2 -10.7 -9.9 -2.3 1.2 -0.4 -0.7 1.4 0.0
0.4 -9.5 -9.4 -2.3 1.2 -0.3 -0.2 1.3 0.2
0.6 -9.5 -9.2 -2.3 1.2 -0.3 -0.9 1.3 0.7
0.8 -10.5 -9.7 -2.3 1.2 -0.3 -0.8 1.0 0.4
"""

etot_ref = {
    "Time": np.array([0.0, 0.2, 0.4, 0.6, 0.8]),
    "ETOT": np.array([-9.7, -10.7, -9.5, -9.5, -10.5]),
    "Madelung": np.array([-9.8, -9.9, -9.4, -9.2, -9.7]),
    "Eigenvalues-Core": np.array([-2.3, -2.3, -2.3, -2.3, -2.3]),
    "Eigenvalues-Valence": np.array([1.2, 1.2, 1.2, 1.2, 1.2]),
    "Exchange": np.array([-0.3, -0.4, -0.3, -0.3, -0.3]),
    "Correlation": np.array([-0.2, -0.7, -0.2, -0.9, -0.8]),
    "XC-potential": np.array([1.3, 1.4, 1.3, 1.3, 1.0]),
    "Coulomb pot. energy": np.array([0.4, 0.0, 0.2, 0.7, 0.4]),
}


def test_parse_etot(tmp_path) -> None:
    file_path = tmp_path / "TOTENERGY_RTTDDFT.OUT"
    file_path.write_text(etot_str)
    etot = parse_etot(file_path.as_posix())
    keys = etot.keys()
    assert keys == etot_ref.keys()
    assert all([np.allclose(etot[key], etot_ref[key]) for key in keys])


forces_str = """ 0.0 0.54  0.82 -0.96
 0.2 0.17 -0.29 -0.52
 0.4 0.15 -0.34 -0.60
 0.6 0.69 -0.39 -0.15
"""


forces_ref = {
    "Time": np.array([0.0, 0.2, 0.4, 0.6]),
    "Fx": np.array([0.54, 0.17, 0.15, 0.69]),
    "Fy": np.array([0.82, -0.29, -0.34, -0.39]),
    "Fz": np.array([-0.96, -0.52, -0.60, -0.15]),
}


def test_parse_force(tmp_path) -> None:
    for name in ["FCR", "FEXT", "FHF", "FVAL"]:
        fname = name + "_0001.OUT"
        file_path = tmp_path / fname
        file_path.write_text(forces_str)
        forces = parse_force(file_path.as_posix())
        keys = forces.keys()
        assert keys == forces_ref.keys()
        assert all([np.allclose(forces[key], forces_ref[key]) for key in keys])


current_str = """ 0.0 0.0 0.0 0.0
0.2  1.6  8.6 -4.3
0.4  0.0 -4.8 -7.1
0.6 -7.4  3.1 -2.2
"""

current_ref = {
    "Time": np.array([0.0, 0.2, 0.4, 0.6]),
    "Jx": np.array([0.0, 1.6, 0.0, -7.4]),
    "Jy": np.array([0.0, 8.6, -4.8, 3.1]),
    "Jz": np.array([0.0, -4.3, -7.1, -2.2]),
}


def test_parse_jind(tmp_path) -> None:
    file_path = tmp_path / "CURRENT.OUT"
    file_path.write_text(current_str)
    current = parse_jind(file_path.as_posix())
    keys = current.keys()
    assert keys == current_ref.keys()
    assert all([np.allclose(current[key], current_ref[key]) for key in keys])


nexc_str = """ Time   N.Elec.GS     N.XS    Sum
0.0 10.0 0.0 10.0
0.2  9.9 0.1 10.0
0.4  9.7 0.3 10.0
0.6  6.5 3.5 10.0
"""

nexc_ref = {
    "Time": np.array([0.0, 0.2, 0.4, 0.6]),
    "number_electrons_GroundState": np.array([10.0, 9.9, 9.7, 6.5]),
    "number_electrons_ExcitedState": np.array([0.0, 0.1, 0.3, 3.5]),
    "sum": np.array([10.0, 10.0, 10.0, 10.0]),
}


def test_parse_nexc(tmp_path) -> None:
    file_path = tmp_path / "N_EXCITATIONS.OUT"
    file_path.write_text(nexc_str)
    nexc = parse_nexc(file_path.as_posix())
    keys = nexc.keys()
    assert keys == nexc_ref.keys()
    assert all([np.allclose(nexc[key], nexc_ref[key]) for key in keys])


occupations_str = """ik =       1
      1     1.97843104
      2     1.10327541
      3     0.57293646
      4     0.31164780
      5     0.21001308
      6     0.00093199

ik =       2
      1     1.99998297
      2     1.99990887
      3     0.00131209
      4     0.00146791
      5     0.00079852
      6     0.00095539

"""

occ_ref = {
    "ik": [1, 2],
    "occupations": [
        np.array([1.97843104, 1.10327541, 0.57293646, 0.31164780, 0.21001308, 0.00093199]),
        np.array([1.99998297, 1.99990887, 0.00131209, 0.00146791, 0.00079852, 0.00095539]),
    ],
}


def test_parse_occupations(tmp_path) -> None:
    file_path = tmp_path / "OCCSV_TXT_10.OUT"
    file_path.write_text(occupations_str)
    occ = parse_occupations(file_path.as_posix())
    assert occ["ik"] == occ_ref["ik"]
    np.testing.assert_allclose(occ["occupations"], occ_ref["occupations"])


proj_file_str_square_matrices = """ ik:          1
   1.00000   0.00000
   0.00000   1.00000
 ik:          2
   1.00000   0.00000
   0.00000   1.00000
"""

reference_parsed_proj_square_matrices = {
    "ik": [1, 2],
    "projection": [np.array([[1.0, 0.0], [0.0, 1.0]]), np.array([[1.0, 0.0], [0.0, 1.0]])],
}

proj_file_str_rectangular_matrices = """ ik:          1
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000
 ik:          2
   0.60000   0.80000   0.00000
   0.00000   0.00000   1.00000
"""

reference_parsed_proj_rectangular_matrices = {
    "ik": [1, 2],
    "projection": [np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]), np.array([[0.6, 0.8, 0.0], [0.0, 0.0, 1.0]])],
}


@pytest.mark.parametrize(
    ["proj_file_str", "reference_parsed_dict"],
    [
        (proj_file_str_square_matrices, reference_parsed_proj_square_matrices),
        (proj_file_str_rectangular_matrices, reference_parsed_proj_rectangular_matrices),
    ],
)
def test_parse_proj_screenshots(proj_file_str, reference_parsed_dict, tmp_path) -> None:
    proj_file_path = tmp_path / "PROJECTION_COEFFS_0.OUT"
    proj_file_path.write_text(proj_file_str)
    proj_out = parse_proj_screenshots(proj_file_path.as_posix())
    is_equal = proj_out["ik"] == reference_parsed_dict["ik"]
    key = "projection"
    is_equal = is_equal and all([np.allclose(x, y) for (x, y) in zip(proj_out[key], reference_parsed_dict[key])])
    assert is_equal


polarization_str = """    0.000      0.000000000000      0.000000000000      0.000000000000
    0.250      0.000042635973      0.000299503212     -0.000231974189
    0.500      0.000147610141      0.000898510863     -0.000695923325
    0.750      0.000169024762      0.001497518242     -0.001159871764
    1.000      0.000183708164      0.002096473012     -0.001623779322
    1.250      0.000347417937      0.002695405154     -0.002087664441
    1.500      0.000396776702      0.003294365230     -0.002551571230
"""

polarization_ref = {
    "Time": np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]),
    "Px": np.array(
        [0.0, 0.000042635973, 0.000147610141, 0.000169024762, 0.000183708164, 0.000347417937, 0.000396776702]
    ),
    "Py": np.array(
        [0.0, 0.000299503212, 0.000898510863, 0.001497518242, 0.002096473012, 0.002695405154, 0.003294365230]
    ),
    "Pz": np.array(
        [0.0, -0.000231974189, -0.000695923325, -0.001159871764, -0.001623779322, -0.002087664441, -0.002551571230]
    ),
}


def test_parse_rttddft_polarization(tmp_path) -> None:
    pol_file_path = tmp_path / "POLARIZATION_RTTDDFT.OUT"
    pol_file_path.write_text(polarization_str)
    pol_out = parse_rttddft_polarization(pol_file_path.as_posix())
    for key in ["Time", "Px", "Py", "Pz"]:
        np.testing.assert_allclose(pol_out[key], polarization_ref[key])
