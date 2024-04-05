"""
Test bse parsers.
The infoxs reference is boiled down to the necessary parts for a test.
Execute tests from exciting_tools directory:
pytest --capture=tee-sys
"""

import numpy as np
import pytest

from excitingtools.exciting_dict_parsers.bse_parser import (
    parse_fastBSE_absorption_spectrum_out,
    parse_fastBSE_exciton_energies_out,
    parse_infoxs_out,
)
from excitingtools.utils.test_utils import MockFile

infoxs_file_str_success = """================================================================================
| EXCITING NITROGEN-14 started for task xsgeneigvec (301)                      =
| version hash id: 1775bff4453c84689fb848894a9224f155377cfc                    =
|                                                                              =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================

 Timings: 
     Date (DD-MM-YYYY)      : 10-12-2020
     Time (hh:mm:ss)        : 20:05:23
     CPU time               : 14.57 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 15 s )
     wall time              : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load               : 684.78 %
     CPU time  (cumulative) : 111.58 sec; 0.03 hrs; ( 0 d, 00 h, 01 m, 52 s )
     wall time (cumulative) : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load  (cumulative) : 684.78 %
 
================================================================================
= EXCITING NITROGEN-14 stopped for task    301                                 =
================================================================================
 
 
================================================================================
| EXCITING NITROGEN-14 started for task writepmatxs (320)                      =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
  Timings: 
     Date (DD-MM-YYYY)      : 10-12-2020
     Time (hh:mm:ss)        : 20:05:23
     CPU time               : 14.57 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 15 s )
     wall time              : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load               : 684.78 %
     CPU time  (cumulative) : 111.58 sec; 0.03 hrs; ( 0 d, 00 h, 01 m, 52 s )
     wall time (cumulative) : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load  (cumulative) : 684.78 %
 
================================================================================
= EXCITING NITROGEN-14 stopped for task    320                                 =
================================================================================
 
 
================================================================================
| EXCITING NITROGEN-14 started for task scrgeneigvec (401)                     =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================

  Timings: 
     Date (DD-MM-YYYY)      : 10-12-2020
     Time (hh:mm:ss)        : 20:05:23
     CPU time               : 14.57 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 15 s )
     wall time              : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load               : 684.78 %
     CPU time  (cumulative) : 111.58 sec; 0.03 hrs; ( 0 d, 00 h, 01 m, 52 s )
     wall time (cumulative) : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load  (cumulative) : 684.78 %
 
================================================================================
= EXCITING NITROGEN-14 stopped for task    401                                 =
================================================================================
 
================================================================================
| EXCITING NITROGEN-14 started for task scrwritepmat (420)                     =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
 
  Timings: 
     Date (DD-MM-YYYY)      : 10-12-2020
     Time (hh:mm:ss)        : 20:05:23
     CPU time               : 14.57 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 15 s )
     wall time              : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load               : 684.78 %
     CPU time  (cumulative) : 111.58 sec; 0.03 hrs; ( 0 d, 00 h, 01 m, 52 s )
     wall time (cumulative) : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load  (cumulative) : 684.78 %

================================================================================
= EXCITING NITROGEN-14 stopped for task    420                                 =
================================================================================
 
================================================================================
| EXCITING NITROGEN-14 started for task bse (445)                              =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================

  Timings: 
     Date (DD-MM-YYYY)      : 10-12-2020
     Time (hh:mm:ss)        : 20:05:23
     CPU time               : 14.57 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 15 s )
     wall time              : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load               : 684.78 %
     CPU time  (cumulative) : 111.58 sec; 0.03 hrs; ( 0 d, 00 h, 01 m, 52 s )
     wall time (cumulative) : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load  (cumulative) : 684.78 %

================================================================================
= EXCITING NITROGEN-14 stopped for task    445                                 =
================================================================================
"""

infoxs_file_str_fail = """================================================================================
| EXCITING NITROGEN-14 started for task xsgeneigvec ( 301)                      =
| version hash id: 1775bff4453c84689fb848894a9224f155377cfc                    =
|                                                                              =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================

  Timings: 
     Date (DD-MM-YYYY)      : 10-12-2020
     Time (hh:mm:ss)        : 20:05:23
     CPU time               : 14.57 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 15 s )
     wall time              : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load               : 684.78 %
     CPU time  (cumulative) : 111.58 sec; 0.03 hrs; ( 0 d, 00 h, 01 m, 52 s )
     wall time (cumulative) : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load  (cumulative) : 684.78 %

================================================================================
= EXCITING NITROGEN-14 stopped for task    301                                 =
================================================================================

================================================================================
| EXCITING NITROGEN-14 started for task writepmatxs ( 320)                     =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================

================================================================================
| EXCITING NITROGEN-14 started for task xsgeneigvec (301)                      =
| version hash id: 1775bff4453c84689fb848894a9224f155377cfc                    =
|                                                                              =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================

  Timings: 
     Date (DD-MM-YYYY)      : 10-12-2020
     Time (hh:mm:ss)        : 20:05:23
     CPU time               : 14.57 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 15 s )
     wall time              : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load               : 684.78 %
     CPU time  (cumulative) : 111.58 sec; 0.03 hrs; ( 0 d, 00 h, 01 m, 52 s )
     wall time (cumulative) : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load  (cumulative) : 684.78 %

================================================================================
= EXCITING NITROGEN-14 stopped for task    301                                 =
================================================================================


================================================================================
| EXCITING NITROGEN-14 started for task writepmatxs (320)                      =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
 
  Timings: 
     Date (DD-MM-YYYY)      : 10-12-2020
     Time (hh:mm:ss)        : 20:05:23
     CPU time               : 14.57 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 15 s )
     wall time              : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load               : 684.78 %
     CPU time  (cumulative) : 111.58 sec; 0.03 hrs; ( 0 d, 00 h, 01 m, 52 s )
     wall time (cumulative) : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load  (cumulative) : 684.78 %
 
================================================================================
= EXCITING NITROGEN-14 stopped for task    320                                 =
================================================================================
 
 
================================================================================
| EXCITING NITROGEN-14 started for task scrgeneigvec (401)                     =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
 
  Timings: 
     Date (DD-MM-YYYY)      : 10-12-2020
     Time (hh:mm:ss)        : 20:05:23
     CPU time               : 14.57 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 15 s )
     wall time              : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load               : 684.78 %
     CPU time  (cumulative) : 111.58 sec; 0.03 hrs; ( 0 d, 00 h, 01 m, 52 s )
     wall time (cumulative) : 2.13 sec; 0.00 hrs; ( 0 d, 00 h, 00 m, 02 s )
     CPU load  (cumulative) : 684.78 %

================================================================================
= EXCITING NITROGEN-14 stopped for task    401                                 =
================================================================================
 
 
================================================================================
| EXCITING NITROGEN-14 started for task scrwritepmat (420)                     =
| Date (DD-MM-YYYY) : 10-12-2020                                               =
================================================================================
"""


reference_parsed_infoxs_file_success = {
    "tasks": [
        {"name": "xsgeneigvec", "number": 301, "finished": True},
        {"name": "writepmatxs", "number": 320, "finished": True},
        {"name": "scrgeneigvec", "number": 401, "finished": True},
        {"name": "scrwritepmat", "number": 420, "finished": True},
        {"name": "bse", "number": 445, "finished": True},
    ],
    "success": True,
    "last_finished_task": "bse",
}

reference_parsed_infoxs_file_fail = {
    "tasks": [
        {"name": "xsgeneigvec", "number": 301, "finished": True},
        {"name": "writepmatxs", "number": 320, "finished": False},
        {"name": "xsgeneigvec", "number": 301, "finished": True},
        {"name": "writepmatxs", "number": 320, "finished": True},
        {"name": "scrgeneigvec", "number": 401, "finished": True},
        {"name": "scrwritepmat", "number": 420, "finished": False},
    ],
    "success": False,
    "last_finished_task": "scrgeneigvec",
}


@pytest.mark.parametrize(
    ["infoxs_file_str", "reference_parsed_dict"],
    [
        (infoxs_file_str_success, reference_parsed_infoxs_file_success),
        (infoxs_file_str_fail, reference_parsed_infoxs_file_fail),
    ],
)
def test_parse_info_xs_out(infoxs_file_str, reference_parsed_dict, tmp_path):
    infoxs_file_path = tmp_path / "INFOXS.OUT"
    infoxs_file_path.write_text(infoxs_file_str)
    info_xs_out = parse_infoxs_out(infoxs_file_path.as_posix())
    assert info_xs_out == reference_parsed_dict


reference_parsed_infoxs_file_times_success = {
    "tasks": [
        {
            "name": "xsgeneigvec",
            "number": 301,
            "finished": True,
            "cpu_time": 14.57,
            "wall_time": 2.13,
            "cpu_time_cum": 111.58,
            "wall_time_cum": 2.13,
        },
        {
            "name": "writepmatxs",
            "number": 320,
            "finished": True,
            "cpu_time": 14.57,
            "wall_time": 2.13,
            "cpu_time_cum": 111.58,
            "wall_time_cum": 2.13,
        },
        {
            "name": "scrgeneigvec",
            "number": 401,
            "finished": True,
            "cpu_time": 14.57,
            "wall_time": 2.13,
            "cpu_time_cum": 111.58,
            "wall_time_cum": 2.13,
        },
        {
            "name": "scrwritepmat",
            "number": 420,
            "finished": True,
            "cpu_time": 14.57,
            "wall_time": 2.13,
            "cpu_time_cum": 111.58,
            "wall_time_cum": 2.13,
        },
        {
            "name": "bse",
            "number": 445,
            "finished": True,
            "cpu_time": 14.57,
            "wall_time": 2.13,
            "cpu_time_cum": 111.58,
            "wall_time_cum": 2.13,
        },
    ],
    "success": True,
    "last_finished_task": "bse",
}

reference_parsed_infoxs_file_times_fail = {
    "tasks": [
        {
            "name": "xsgeneigvec",
            "number": 301,
            "finished": True,
            "cpu_time": 14.57,
            "wall_time": 2.13,
            "cpu_time_cum": 111.58,
            "wall_time_cum": 2.13,
        },
        {"name": "writepmatxs", "number": 320, "finished": False},
        {
            "name": "xsgeneigvec",
            "number": 301,
            "finished": True,
            "cpu_time": 14.57,
            "wall_time": 2.13,
            "cpu_time_cum": 111.58,
            "wall_time_cum": 2.13,
        },
        {
            "name": "writepmatxs",
            "number": 320,
            "finished": True,
            "cpu_time": 14.57,
            "wall_time": 2.13,
            "cpu_time_cum": 111.58,
            "wall_time_cum": 2.13,
        },
        {
            "name": "scrgeneigvec",
            "number": 401,
            "finished": True,
            "cpu_time": 14.57,
            "wall_time": 2.13,
            "cpu_time_cum": 111.58,
            "wall_time_cum": 2.13,
        },
        {"name": "scrwritepmat", "number": 420, "finished": False},
    ],
    "success": False,
    "last_finished_task": "scrgeneigvec",
}


@pytest.mark.parametrize(
    ["infoxs_file_str", "reference_parsed_dict"],
    [
        (infoxs_file_str_success, reference_parsed_infoxs_file_times_success),
        (infoxs_file_str_fail, reference_parsed_infoxs_file_times_fail),
    ],
)
def test_parse_info_xs_out_timing(infoxs_file_str, reference_parsed_dict, tmp_path):
    infoxs_file_path = tmp_path / "INFOXS.OUT"
    infoxs_file_path.write_text(infoxs_file_str)
    info_xs_out = parse_infoxs_out(infoxs_file_path.as_posix(), parse_timing=True)
    assert info_xs_out == reference_parsed_dict


@pytest.fixture
def fastBSE_absorption_spectrum_out_mock(tmp_path):
    string_contents = """# fastBSE imaginary macroscopic dielectric function
# 
# Energy unit:  0.3674932539796232E-01 Hartree
# Broadening:   0.1360569193000000E+00 energy unit
#
#                 omega                    oc11                    oc22                    oc33
+0.0000000000000000E+00 +0.0000000000000000E+00 +0.0000000000000000E+00 +0.0000000000000000E+00
+0.6046974191111111E+01 +0.2462935121365850E+00 +0.2462935121365850E+00 +0.2462935121365850E+00
+0.1209394838222222E+02 +0.2028366299352899E+02 +0.2028366299352899E+02 +0.2028366299352899E+02"""

    file = tmp_path / "fastBSE_absorption_spectrum.out"
    file.write_text(string_contents)
    return MockFile(file, string_contents)


reference_fastBSE_absorption_spectrum = {
    "energy_unit": 0.3674932539796232e-01,
    "broadening": 0.1360569193000000e00,
    "frequency": np.array([+0.0000000000000000e00, +0.6046974191111111e01, +0.1209394838222222e02]),
    "imag_epsilon": np.array(
        [
            [+0.0000000000000000e00, +0.0000000000000000e00, +0.0000000000000000e00],
            [+0.2462935121365850e00, +0.2462935121365850e00, +0.2462935121365850e00],
            [+0.2028366299352899e02, +0.2028366299352899e02, +0.2028366299352899e02],
        ]
    ),
}


def test_parse_fastBSE_absorption_spectrum_out_parser(fastBSE_absorption_spectrum_out_mock):
    fastBSE_absorption_spectrum_out = parse_fastBSE_absorption_spectrum_out(fastBSE_absorption_spectrum_out_mock.file)
    assert np.allclose(
        fastBSE_absorption_spectrum_out["energy_unit"], reference_fastBSE_absorption_spectrum["energy_unit"]
    )
    assert np.allclose(
        fastBSE_absorption_spectrum_out["broadening"], reference_fastBSE_absorption_spectrum["broadening"]
    )
    assert np.allclose(fastBSE_absorption_spectrum_out["frequency"], reference_fastBSE_absorption_spectrum["frequency"])
    assert np.allclose(
        fastBSE_absorption_spectrum_out["imag_epsilon"], reference_fastBSE_absorption_spectrum["imag_epsilon"]
    )


@pytest.fixture
def fastBSE_exciton_energies_out_mock(tmp_path):
    string_contents = """# fastBSE exciton eigen energies
# The three rows correspond to the results of the three Lanczos runs, each for one of the
# directions of <p> as starting point.
# 
# Energy unit:  0.3674932539796232E-01 Hartree
# IP band gap:  0.8205751967708917E+01 energy unit
# 
#            E -> <p_1>              E -> <p_2>              E -> <p_3>
+0.8166524450038512E+01 +0.8166363823856052E+01 +0.8166439398878977E+01
+0.8168631668065551E+01 +0.8198150340644158E+01 +0.8588070401753557E+01
+0.8587998131087492E+01 +0.8588200092979093E+01 +0.9161347631625411E+01"""

    file = tmp_path / "fastBSE_exciton_energies.out"
    file.write_text(string_contents)
    return MockFile(file, string_contents)


reference_fastBSE_exciton_energies = {
    "energy_unit": 0.3674932539796232e-01,
    "ip_band_gap": 0.8205751967708917e01,
    "exciton_energies": np.array(
        [
            [+0.8166524450038512e01, +0.8166363823856052e01, +0.8166439398878977e01],
            [+0.8168631668065551e01, +0.8198150340644158e01, +0.8588070401753557e01],
            [+0.8587998131087492e01, +0.8588200092979093e01, +0.9161347631625411e01],
        ]
    ),
}


def test_parse_fastBSE_exciton_energies_out_parser(fastBSE_exciton_energies_out_mock):
    fastBSE_exciton_energies_out = parse_fastBSE_exciton_energies_out(fastBSE_exciton_energies_out_mock.file)
    assert np.allclose(fastBSE_exciton_energies_out["energy_unit"], reference_fastBSE_exciton_energies["energy_unit"])
    assert np.allclose(fastBSE_exciton_energies_out["ip_band_gap"], reference_fastBSE_exciton_energies["ip_band_gap"])
    assert np.allclose(
        fastBSE_exciton_energies_out["exciton_energies"], reference_fastBSE_exciton_energies["exciton_energies"]
    )
