"""
Module for functions that determine job dependencies and batching
"""
from typing import List, Dict


def compute_dependency_trees(dependencies: dict) -> Dict[str, List[str]]:
    """
    For each job, determine all prior jobs that need to run before it can be executed.

    Each key of 'dependencies' defines a job, and each value represents its
    immediate (or explicit) dependencies. This routine recursively looks through the
    dependencies of dependencies, such that a list of ALL direct and indirect
    (or explicit and implicit) dependencies is returned, for each job.

    See test_functions/test_dependencies for clear examples.

    If a circular dependency is found, an exception is immediately raised as it implies
    an error in the configuration file that defines the test cases.

    :params dict dependencies: Dependencies for each job, where key = job name
    and value = list of dependencies.
    :return Dict[str, List[str]] all_dependencies: Jobs and their total dependencies.
    """

    # Recursively find all dependencies for a given key of `dependencies`
    def get_all_dependencies(dependencies_of_key: list, combined_deps: list) -> List[str]:
        if dependencies_of_key == ['']:
            return combined_deps

        for d in dependencies_of_key:
            combined_deps.append(d)
            get_all_dependencies(dependencies[d], combined_deps)

        return combined_deps

    # Dict of all dependencies (explicit and implicit) per key
    all_dependencies = {}

    for key, explicit_dependencies_of_key in dependencies.items():
        all_dependencies_of_key: list = []
        try:
            all_dependencies[key] = get_all_dependencies(explicit_dependencies_of_key,
                                                         all_dependencies_of_key)
        except RecursionError:
            raise RecursionError(f'Test {key} has circular dependencies')

    return all_dependencies


def schedule_tests(dependencies: dict) -> dict:
    """
    Schedule tests such that one runs in batches, defined by jobs with:
    0, 1, 2, ... N dependencies.
    This guarantees all dependencies are run prior to a job with N dependencies running.

    Batch keys with empty values are excluded prior to returning.
    For example, batch[i] = [] , for any i, will be removed.

    :param dict dependencies: Jobs and associated dependencies.
    :return dict batch: Jobs ordered according to order of execution, where key = execution order
    and value = tests associated with the execution order.
    """
    max_dependencies = max([len(deps) for deps in dependencies.values()])
    batch = {i: [] for i in range(0, max_dependencies + 1)}

    for test_name, deps in dependencies.items():
        batch[len(deps)].append(test_name)

    return {i: deps for i, deps in batch.items() if deps}
