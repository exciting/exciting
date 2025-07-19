"""Test of get_hall_number_from_symmetry."""

from __future__ import annotations

import pytest
from spglib import get_hall_number_from_symmetry


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_get_hall_number_from_symmetry(crystal_data_dataset):
    """Test get_hall_number_from_symmetry."""
    dataset = crystal_data_dataset["dataset"]
    hall_number = get_hall_number_from_symmetry(
        dataset.rotations,
        dataset.translations,
        symprec=crystal_data_dataset["symprec"],
    )
    assert hall_number == dataset.hall_number
