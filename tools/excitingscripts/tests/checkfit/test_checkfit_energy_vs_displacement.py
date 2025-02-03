import json
import os
import warnings

import numpy as np
from excitingscripts.checkfit.checkfit import quantity_specific_checkfit

try:
    # Check if np.RankWarning exists directly
    numpy_RankWarning = np.RankWarning
except AttributeError:
    # If not, fall back to np.exceptions.RankWarning
    numpy_RankWarning = np.exceptions.RankWarning

def compare_arrays_with_none(array, ref_array):
    # Ensure both lists have the same length
    assert len(array) == len(ref_array), "List lengths do not match"

    # Filter out None values for both lists and compare the rest
    for arr, ref_arr in zip(array, ref_array):
        arr_filtered = [a for a in arr if a is not None]
        ref_arr_filtered = [r for r in ref_arr if r is not None]

        # Compare filtered lists
        assert np.allclose(np.array(arr_filtered), np.array(ref_arr_filtered))

def test_checkfit_energy_vs_displacement(phonon_results_mock, info_diamond_phonon_mock, tmp_path):
    phonon_results_file = tmp_path / "phonon_results.json"
    phonon_results_file.write_text(phonon_results_mock.string)
    info_diamond_phonon_file = tmp_path / "INFO-diamond-phonon.json"
    info_diamond_phonon_file.write_text(info_diamond_phonon_mock.string)

    os.chdir(tmp_path)
    warnings.simplefilter("error", numpy_RankWarning)
    check_fit_func = quantity_specific_checkfit("energy", 2 / 3, 2)
    check_fit_func(0.025, 2, 12.01)

    checkfit_results = tmp_path / "checkfit_energy_results.json"
    with open(checkfit_results) as fid:
        file_content = json.load(fid)

    ref_order_of_derivative = 2
    assert file_content["order_of_derivative"] == ref_order_of_derivative

    ref_max_displacement_values = [0.025, 0.02, 0.015, 0.01, 0.005]
    ref_frequencies = [[1602.114306187767, 1602.115393552411, 1584.1696938890245, 1584.1696618360838, 1584.04486817542,
                        1584.0448720964116],
                       [1595.9906790722566, 1595.9920142789936, 1584.0987032841483, 1584.0986795448619,
                        1583.9958213776597, 1583.9958200976864],
                       [1591.0415284637754, 1591.0432682904334, 1584.0256415215833, 1584.0256222867367,
                        1583.9829277587307, None],
                       [1587.2747592698036, 1587.277291091997, 1583.9872707841052, None, None, None],
                       [1584.7255547086195, None, None, None, None, None]]

    max_displacement_values = []
    frequencies = []
    for data in file_content["fits"]:
        max_displacement_values.append(data["max_displacement"])
        frequencies.append(data["frequencies"])

    assert np.allclose(np.array(max_displacement_values), np.array(ref_max_displacement_values))
    compare_arrays_with_none(frequencies, ref_frequencies)
