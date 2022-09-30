from excitingtools.dataclasses.band_structure import BandData
import pytest
import numpy as np


@pytest.fixture
def band_data():
    """ Initialise an instance of BandData, and check attributes are set correctly.
    Reference data taken from 'bandstructure.xml' and 'bandstructure.dat' files for silver, containing only two bands
    and only 6 k-sampling points per band.
    """

    ref_bands = np.array([[-3.37071333, 0.59519479],
                          [-3.37071074, 0.59537998],
                          [-3.37070519, 0.59575556],
                          [-3.3706986, 0.59616303],
                          [-3.3706822, 0.59664939],
                          [-3.37066123, 0.59639211]])

    ref_k_points = np.array([[1.000000, 0.000000, 0.000000],
                             [0.988281, 0.011719, 0.000000],
                             [0.976562, 0.023438, 0.000000],
                             [0.964844, 0.035156, 0.000000],
                             [0.953125, 0.046875, 0.000000],
                             [0.941406, 0.058594, 0.000000]])

    ref_e_fermi = 0.0

    ref_flattened_k_points = np.array([0., 0.02697635, 0.05395270, 0.08092905, 0.10790540, 0.13488176])

    ref_vertices = [{'distance': 0.0, 'label': 'Gamma', 'coord': [0.0, 0.0, 0.0]},
                    {'distance': 0.8632432750, 'label': 'K', 'coord': [0.625, 0.375, 0.0]},
                    {'distance': 1.150991033, 'label': 'X', 'coord': [0.5, 0.5, 0.0]},
                    {'distance': 1.964864598, 'label': 'G', 'coord': [0.0, 0.0, 0.0]},
                    {'distance': 2.669699781, 'label': 'L', 'coord': [0.5, 0.0, 0.0]}]

    band_data = BandData(ref_bands, ref_k_points, ref_e_fermi, ref_flattened_k_points, ref_vertices)

    assert band_data.n_k_points == band_data.bands.shape[0], (
        "First dim of bands array equals the number of k-sampling points in the band structure")
    assert band_data.n_k_points == 6, "sampling points per band"
    assert band_data.n_bands == 2, "band_structure_xml contains two bands"
    assert np.allclose(band_data.k_points, ref_k_points, atol=1.e-8)
    assert np.allclose(band_data.bands, ref_bands, atol=1.e-8)
    assert band_data.vertices == ref_vertices

    return band_data


def test_xticks_and_labels(band_data):
    flattened_high_sym_points_ref = [0.000000000, 0.86324327, 1.15099103, 1.9648646, 2.66969978]
    unicode_gamma = '\u0393'
    labels_ref = [unicode_gamma, "K", "X", unicode_gamma, "L"]

    flattened_high_sym_points, labels = band_data.band_path()

    assert np.allclose(flattened_high_sym_points, flattened_high_sym_points_ref)
    assert labels == labels_ref


def test_get_band_edges(band_data):
    i_vbm, i_cbm = band_data.get_band_edges()
    assert i_vbm == 0
    assert i_cbm == 1

    band_data.e_fermi = 1.0
    with pytest.raises(ValueError, match="Fermi level 1.0 larger than highest band energy 0.59664939"):
        band_data.get_band_edges(), "ValueError is returned in case of erroneously large Fermi level"


def test_get_valence_band_maximum(band_data):
    vbm = band_data.get_valence_band_maximum()
    assert vbm == -3.37066123


def test_get_conduction_band_minimum(band_data):
    cbm = band_data.get_conduction_band_minimum()
    assert cbm == 0.59519479


def test_get_band_gap(band_data):
    k_valence = [1.000000, 0.000000, 0.000000]
    k_conduction = [0.941406, 0.058594, 0.000000]
    band_gap = band_data.get_band_gap(k_valence, k_conduction)
    assert band_gap == 3.96710544


def test_get_fundamental_band_gap(band_data):
    fundamental_gap = band_data.get_fundamental_band_gap()
    assert fundamental_gap == 3.96585602

    band_data.bands[0, 0] = 1
    fundamental_gap_zero = band_data.get_fundamental_band_gap()
    assert fundamental_gap_zero == 0.0
