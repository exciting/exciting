""" Run references for all test cases.

There are very few use cases where this should be done!
Note, this has not been tested following refactoring into its own module.
"""
import os

from ..exciting_settings.constants import settings as exciting_settings, Defaults
from ..reference.reference import run_single_reference
from ..runner.execute import set_job_environment
from ..runner.profile import build_type_enum_to_str
from ..runner.set_tests import get_test_cases_from_config


def run_references(settings: Defaults):
    """Run all test suite cases, such that all reference outputs are regenerated.
    There are very few use cases where this should be done

    :param settings: Default settings for environment variables, paths and run-time parameters.
    """
    executable = os.path.join(settings.exe_dir, build_type_enum_to_str[settings.exe_ref])
    my_env = set_job_environment({'omp': '1', 'mkl_threads': '1'})

    tests_in_config = get_test_cases_from_config(settings)
    print(f'Running tests defined in {settings.tests_file}')
    print(f"All references will be run with {settings.exe_ref}")

    for test in list(tests_in_config):
        reference_dir = os.path.join(test, settings.ref_dir)
        run_single_reference(reference_dir, executable, settings, my_env)


if __name__ == "__main__":
    while True:
        answer = input("Are you sure that you want to rerun references? "
                       "This will replace ALL existing references with references "
                       "created by the current state of the code and should only "
                       "be done if all tests succeed with the current version of the code.(yes/no)\n")
        if answer.lower() in ['yes', 'y']:
            break
        elif answer.lower() in ['no', 'n']:
            exit()
    run_references(exciting_settings)
