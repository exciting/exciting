"""
 Failing and hanging tests.

 Each skipped test should have a TODO(Name) assigned to it here, which
 is also logged on the gitlab issue tracker.

 Each list entry should be a dictionary.
"""
import sys

from .utils import Compiler, Build_type, CompilerBuild, get_compiler_type, build_type_str_to_enum

# Flakey tests Issue #54
# Update 16th March 2020. Remove the hybrid tests from the suite until more
# robust testing is introduced and someone investigates them 
#
# Tests that don't always pass the CI for any discernable reason
# are considered flakey rather than failing.
#
# This currently corresponds to the hybrids:
#  - groundstate-HYB_HSE-Si
#  - groundstate-HYB_PBE0-Si
#
# When ctest is introduced, this can addressed with the -repeat until-pass:2 option.
# Additionly, explicit assertions and tolerances will make the hybrid tests more robust.
# See MR !38 Testing/pytest.
#
# Also see Issue 51. Review tolerances for all tests migrated from the old test suite.


failing_tests = [
    # TODO(Alex) Issue #35. MPI GW calculations do not produce EPS00_GW.OUT
    # Also appear to get a failure for Intel 2019 serial (only occurs on the CI)
    {'name': 'GW/LDA_PW-gw-Si',
     'comment': 'MPI GW calculations do not produce EPS00_GW.OUT',
     'tags': [CompilerBuild(Compiler.intel, Build_type.serial),
              CompilerBuild(Compiler.all, Build_type.mpiandsmp)]
     },

    # TODO(Sebastian) Issue 40
    {'name': 'properties/LDA_PZ-wannier-SiC',
     'comment': '',
     'tags': [CompilerBuild(Compiler.gcc, Build_type.all)]
     },

    # TODO(Sven) Issue #39
    # Andris comment: I have tried several runs with a different number of threads and MPI ranks,
    # but the calculations always converged in 19-20 iterations and always to the same energy with the threshold.
    # - Looks like issue with comparing integers in the test suite, rather than a test failure
    # - One just needs to set a tolerance for integers
    # {'name': 'groundstate/LDA_PW-noncollinear-Fe',
    #  'comment': 'Number of SCF iterations differs to reference',
    #  'tags': [CompilerBuild(Compiler.all, Build_type.mpiandsmp)]
    #  },

    # TODO ADD ISSUE
    {'name': 'properties/PBE-properties-Si',
     'comment': 'Epsilon 11 and 33 do not agree with serial reference values',
     'tags': [CompilerBuild(Compiler.all, Build_type.mpiandsmp)]
     },

    # TODO(Cecilia) Issue 54 
    {'name': 'hybrids/HSE-Si',
     'comment': 'Test is flakey when run in the CI with Intel parallel build',
     'tags': [CompilerBuild(Compiler.intel, Build_type.mpiandsmp)]
     },

    # TODO(Cecilia) Issue 54 
    {'name': 'hybrids/PBE0-Si',
     'comment': 'Test is flakey when run in the CI with Intel parallel build',
     'tags': [CompilerBuild(Compiler.intel, Build_type.mpiandsmp)]
     },

    # TODO(Maria) Issue 55
    {'name': 'properties/LDA_PW-transport-Si',
     'comment': 'Test is flakey when run in the CI with GCC builds: Test outputs are not written',
     'tags': [CompilerBuild(Compiler.gcc, Build_type.all)]
     },

    # TODO(Alex) Issue 101. Test still passes with Intel MPIandSMP. Also see Issue #75
    {'name': 'groundstate/LDA_PW-collinear-Fe',
     'comment': 'Most energies differ to reference by ~1.e-7. '
                'scl%Sum of eigenvalues by ~ 1.e-6. '
                'DOS at Fermi differs by 5.7e-04. ',
     'tags': [CompilerBuild(Compiler.gcc, Build_type.all)]
     },
    {'name': 'groundstate/LDA_PW-collinear-Fe',
     'comment': 'Most energies differ to reference by ~1.e-7. '
                'scl%Sum of eigenvalues by ~ 1.e-6. '
                'DOS at Fermi differs by 5.7e-04. ',
     'tags': [CompilerBuild(Compiler.intel, Build_type.all)]
     }
]

# These always have to be skipped else the suite hangs 
hanging_tests = [
    # TODO(Alex) Issue #36 chargedensityplot hangs when running with np > 1 cores
    # Should remove chargedensityplot from this and add it as a property
    {'name': 'groundstate/PBE-Al',
     'comment': 'chargedensityplot hangs at an allgatherv call',
     'tags': [CompilerBuild(Compiler.all, Build_type.mpiandsmp)]}
]


def set_skipped_tests(executable: str, incl_failing_tests: bool) -> list:
    """
    Generate full list of tests to skip from the test suite

    :param str executable:
    :param bool incl_failing_tests: bool. If true, include failing tests in the suite,
    hence exclude from the list of skipped tests

    :return list tests_to_skip: list of tests to skip, where each entry is a dict
    """

    compiler = get_compiler_type()

    if compiler is None:
        print('Compiler could not be determined from build/make.inc\n'
              'Excluding all failing tests from the test suite')
        compiler = Compiler.all

    build_type = build_type_str_to_enum[executable]
    tests_to_skip = []

    # Always include hanging tests in tests to skip, else the test suite
    # will hang
    for test in hanging_tests:
        for tag in test['tags']:
            if (tag.compiler == compiler or tag.compiler == Compiler.all) \
                    and (tag.build == build_type or tag.build == Build_type.all):
                tests_to_skip.append(test)

    if incl_failing_tests:
        return tests_to_skip

    for test in failing_tests:
        for tag in test['tags']:
            if (tag.compiler == compiler or tag.compiler == Compiler.all) \
                    and (tag.build == build_type or tag.build == Build_type.all):
                tests_to_skip.append(test)

    return tests_to_skip
