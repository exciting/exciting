"""Test of spacegroup_type_from_symmetry."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
from spglib import (
    get_spacegroup_type,
    get_spacegroup_type_from_symmetry,
)

if TYPE_CHECKING:
    pass


@pytest.mark.parametrize("lattice_None", [True, False])
def test_spacegroup_type_from_symmetry(
    lattice_None: bool,
    crystal_data_dataset,
):
    """Test spacegroup_type_from_symmetry."""
    if lattice_None:
        lattice = None
    else:
        lattice = crystal_data_dataset["crystal_data"].cell[0]

    dataset = crystal_data_dataset["dataset"]
    spgtype = get_spacegroup_type_from_symmetry(
        dataset.rotations,
        dataset.translations,
        lattice=lattice,
        symprec=crystal_data_dataset["symprec"],
    )

    assert spgtype.number == dataset.number
    spgtype_ref = get_spacegroup_type(dataset.hall_number)
    assert spgtype == spgtype_ref
