from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
from spglib import get_symmetry_dataset

if TYPE_CHECKING:
    from conftest import CrystalData


@pytest.mark.benchmark(group="space-group")
def test_get_symmetry_dataset(
    benchmark,
    crystal_data: CrystalData,
):
    """Benchmarking get_symmetry_dataset on all structures under test/data."""
    # TODO: This benchmarks individual crystal data. Is it useful?
    print(f"Benchmark get_symmetry_dataset on {crystal_data.name} structures")

    def _get_symmetry_dataset_for_cells():
        _ = get_symmetry_dataset(crystal_data.cell, symprec=1e-5)

    benchmark.pedantic(_get_symmetry_dataset_for_cells, rounds=4)
