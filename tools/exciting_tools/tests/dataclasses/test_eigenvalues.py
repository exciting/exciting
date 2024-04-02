from typing import List

import numpy as np
import pytest

from excitingtools.dataclasses.data_structs import BandIndices, NumberOfStates, PointIndex
from excitingtools.dataclasses.eigenvalues import EigenValues


@pytest.fixture
def eigenvalues_instance():
    """Initialise an instance of EigenValues, and check attributes are set correctly.
    Reference data taken from `evalqp_oxygen` in test_gw_eigenvalues.py

    G0W0 band structure

    Band index of VBM:  21
    Band index of CBm:  22

    Indirect BandGap (eV):                    5.3790
    at k(VBM) =    0.000   0.500   0.500 ik =     3
       k(CBm) =    0.000   0.000   0.000 ik =     1
    Direct Bandgap at k(VBM) (eV):            5.5451
    Direct Bandgap at k(CBm) (eV):            5.9468
    """
    ref_gw_eigenvalues = np.array(
        [
            [-0.12905, -0.12891, -0.12896, 0.08958, 0.08957, 0.18396],
            [-0.19801, -0.17047, -0.17035, 0.11429, 0.11430, 0.19514],
            [-0.15828, -0.15818, -0.10809, 0.09569, 0.14613, 0.18404],
        ]
    )

    ref_k_points = np.array(
        [[0.000000, 0.000000, 0.000000], [0.000000, 0.000000, 0.500000], [0.000000, 0.500000, 0.500000]]
    )

    ref_weights = [0.125, 0.5, 0.375000]

    ref_occupations = np.array([[2, 2, 2, 0, 0, 0], [2, 2, 2, 0, 0, 0], [2, 2, 2, 0, 0, 0]])

    eigen_values = EigenValues(
        NumberOfStates(19, 24), ref_k_points, [1, 2, 3], ref_gw_eigenvalues, ref_weights, ref_occupations
    )

    assert ref_gw_eigenvalues.shape == (len(eigen_values.k_points), eigen_values.state_range.n_states)
    assert eigen_values.state_range.first_state == 19
    assert eigen_values.state_range.last_state == 24
    assert np.allclose(eigen_values.k_points, ref_k_points)
    assert eigen_values.k_indices == [1, 2, 3]
    assert np.allclose(eigen_values.all_eigenvalues, ref_gw_eigenvalues), "GW column eigenvalues, for all k-points"
    assert eigen_values.weights == ref_weights
    assert np.allclose(eigen_values.occupations, ref_occupations)
    return eigen_values


def test_class_eigenvalues_get_array_index(eigenvalues_instance):
    assert eigenvalues_instance.get_array_index(19) == 0, "Index of first state in array"
    assert eigenvalues_instance.get_array_index(24) == 5, "Index of last state in array"
    assert eigenvalues_instance.state_range.n_states == 6, "Total number of states"


def test_class_eigenvalues_get_k_point(eigenvalues_instance):
    assert np.allclose(eigenvalues_instance.get_k_point(1), [0.000000, 0.000000, 0.000000])
    assert np.allclose(eigenvalues_instance.get_k_point(2), [0.000000, 0.000000, 0.500000])
    assert np.allclose(eigenvalues_instance.get_k_point(3), [0.000000, 0.500000, 0.500000])


def test_class_eigenvalues_get_index(eigenvalues_instance):
    assert eigenvalues_instance.get_index([0.000000, 0.000000, 0.000000]) == 1
    assert eigenvalues_instance.get_index([0.000000, 0.000000, 0.500000]) == 2
    assert eigenvalues_instance.get_index([0.000000, 0.500000, 0.500000]) == 3
    assert eigenvalues_instance.get_index([0.500000, 0.500000, 0.500000]) == EigenValues.NO_MATCH, "No k-point matched"


def test_class_eigenvalues_get_k_points(eigenvalues_instance):
    k_points_and_indices: List[PointIndex] = eigenvalues_instance.get_k_points()
    assert np.allclose(k_points_and_indices[0].point, [0.000000, 0.000000, 0.000000])
    assert np.allclose(k_points_and_indices[1].point, [0.000000, 0.000000, 0.500000])
    assert np.allclose(k_points_and_indices[2].point, [0.000000, 0.500000, 0.500000])
    assert k_points_and_indices[0].index == 1
    assert k_points_and_indices[1].index == 2
    assert k_points_and_indices[2].index == 3


def test_class_eigenvalues_get_eigenvalues(eigenvalues_instance):
    eigenvalues = eigenvalues_instance.get_eigenvalues(k_point=[0.0, 0.5, 0.5])
    assert np.allclose(eigenvalues, [-0.15828, -0.15818, -0.10809, 0.09569, 0.14613, 0.18404])

    assert eigenvalues_instance.get_eigenvalues(k_point=[0.5, 0.5, 0.5]).size == 0, "No k-point, hence eigenvalues"


def test_class_eigenvalues_band_gap(eigenvalues_instance):
    k_valence = [0.000, 0.500, 0.500]
    k_conduction = [0.000, 0.000, 0.000]

    indirect_band_gap = eigenvalues_instance.band_gap(BandIndices(21, 22), k_points=[k_valence, k_conduction])
    direct_band_gap_at_Gamma = eigenvalues_instance.band_gap(BandIndices(21, 22), k_points=[k_conduction, k_conduction])

    assert np.isclose(indirect_band_gap, 0.19767), "Indirect band gap in Ha"
    assert np.isclose(direct_band_gap_at_Gamma, 0.218540887), "Direct band gap at Gamma, in Ha"


def test_class_eigenvalues_get_transition_energy(eigenvalues_instance):
    k_valence = [0.000, 0.500, 0.500]
    k_conduction = [0.000, 0.000, 0.000]

    indirect_gap = eigenvalues_instance.get_transition_energy(k_valence, k_conduction)
    direct_gap_at_Gamma = eigenvalues_instance.get_transition_energy(k_conduction, k_conduction)

    assert np.isclose(indirect_gap, 0.19767), "Transition energy for X→Γ, in Ha"
    assert np.isclose(direct_gap_at_Gamma, 0.21854), "Transition energy for Γ→Γ, in Ha"

    with pytest.raises(ValueError, match="Requested conduction k-point \\[1, 1, 1\\] not present."):
        (
            eigenvalues_instance.get_transition_energy(k_valence, [1, 1, 1]),
            "ValueError is returned if k-point is not matched.",
        )
