"""
Test dependencies module
"""
import pytest
from collections import OrderedDict

from ..src.runner.dependencies import  compute_dependency_trees, schedule_tests


def test_compute_dependency_trees():
    """
    Recursively find all dependencies, per key.
    """

    test_dependencies = {'E': [''],
                         'A': [''],
                         'B': ['E'],
                         'C': ['B', 'A'],  # C depends on A and B, A and B are independent, B depends on E
                         'D': ['C']}  # D depends on C, therefore implicitly depends also on A, B and E

    all_dependencies = compute_dependency_trees(test_dependencies)
    ref_all_dependencies = {'E': [], 'A': [], 'B': ['E'], 'C': ['B', 'E', 'A'], 'D': ['C', 'B', 'E', 'A']}
    assert all_dependencies == ref_all_dependencies


def test_compute_dependency_trees_circular_dependencies():

    test_dependencies = OrderedDict({'E': ['A'],
                                     'A': ['E']}
                                    )

    with pytest.raises(RecursionError) as error_info:
        all_dependencies = compute_dependency_trees(test_dependencies)

    assert error_info.value.args[0] == 'Test E has circular dependencies'


def test_schedule_tests():

    dependencies = {'E': [], 'A': [], 'B': ['E'], 'C': ['B', 'E', 'A'], 'D': ['C', 'B', 'E', 'A']}

    ref_batch = {0: ['E', 'A'], 1: ['B'], 3: ['C'], 4: ['D']}
    batch = schedule_tests(dependencies)

    assert 2 not in batch, "No tests have 2 total dependencies"

    assert batch == ref_batch, \
        "Expect keys with no dependencies in batch[0], keys with 1 dependency in batch[0], " \
        "... keys with N dependencies in batch[N]"
