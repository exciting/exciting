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

def test_checkfit_force_vs_displacement(phonon_results_mock, info_diamond_phonon_mock, tmp_path):
    phonon_results_file = tmp_path / "phonon_results.json"
    phonon_results_file.write_text(phonon_results_mock.string)
    info_diamond_phonon_file = tmp_path / "INFO-diamond-phonon.json"
    info_diamond_phonon_file.write_text(info_diamond_phonon_mock.string)

    os.chdir(tmp_path)
    warnings.simplefilter("error", numpy_RankWarning)
    check_fit_func = quantity_specific_checkfit("force", 2, 1)
    check_fit_func(0.025, 1, 12.01)

    checkfit_results = tmp_path / "checkfit_force_results.json"
    with open(checkfit_results) as fid:
        file_content = json.load(fid)

    ref_order_of_derivative = 1
    assert file_content["order_of_derivative"] == ref_order_of_derivative

    ref_max_displacement_values = [0.025, 0.02, 0.015, 0.01, 0.005]
    ref_frequencies = [[1605.965154684306, 1605.9644361781159, 1580.6004156696795, 1580.6004493734529,
                        1580.622640170721, 1580.622638444908],
                       [1597.4581838372098, 1597.457308333169, 1580.6142617703342, 1580.6142882061579,
                        1580.6163954958577, 1580.6163944798159],
                       [1590.628311307824, 1590.6271881404223, 1580.6160033168942, 1580.6160228968226,
                        1580.6128094418616, 1580.6128090593706],
                       [1585.488274061948, 1585.4867039688875, 1580.6133110956569, 1580.6133242413734, None, None],
                       [1582.0508444148923, 1582.0482299987625, None, None, None, None]]

    max_displacement_values = []
    frequencies = []
    for data in file_content["fits"]:
        max_displacement_values.append(data["max_displacement"])
        frequencies.append(data["frequencies"])

    assert np.allclose(np.array(max_displacement_values), np.array(ref_max_displacement_values))
    compare_arrays_with_none(frequencies, ref_frequencies)
