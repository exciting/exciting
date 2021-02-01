"""
 Failing and hanging tests.

 Each skipped test should have a TODO(Name) assigned to it here, which
 is also logged on the gitlab issue tracker.

 Each list entry should be a dictionary.
"""
import sys 
sys.path.insert(1, 'tools')
from utils import Compiler, Build_type, CompilerBuild, get_compiler_type, build_type_enum_map


failing_tests = [
    # TODO(Alex) Issue #35. MPI GW calculations do not produce EPS00_GW.OUT
    # Also appear to get a failure for Intel 2019 serial (only occurs on the CI)
    {'name':'groundstate-LDA_PW-gw-Si',
     'comment':'MPI GW calculations do not produce EPS00_GW.OUT',
     'tags': [CompilerBuild(Compiler.intel, Build_type.serial), 
              CompilerBuild(Compiler.all, Build_type.mpiandsmp)]
    }, 

    # TODO(Sebastian) Issue 40
    {'name': 'groundstate-LDA_PZ-wannier-SiC',
     'comment':'',
     'tags': [CompilerBuild(Compiler.gcc, Build_type.all)]
    }, 

    # TODO(Ronaldo) Issue #34. 
    # MPI calculations result in numerous files containing numbers that differ from serial calculations
    {'name': 'groundstate-GGA_PBE-tddft-LiF',
     'comment': 'Serial and MPI numbers differ',
     'tags': [CompilerBuild(Compiler.all, Build_type.mpiandsmp)]
    },

    # TODO(Sven) Issue #39 
    # Andris comment: I have tried several runs with a different number of threads and MPI ranks, 
    # but the calculations always converged in 19-20 iterations and always to the same energy with the threshold.
    # - Looks like issue with comparing integers in the test suite, rather than a test failure
    {'name':'groundstate-LDA_PW-noncollinear-Fe',
     'comment': 'Number of SCF iterations differs to reference',
     'tags': [CompilerBuild(Compiler.gcc, Build_type.mpiandsmp)]
     }, 

    # TODO ADD ISSUE
    {'name':'groundstate-GGA_PBE-properties-Si',
     'comment':'Epislon 11 and 33 do not agree with serial reference values',
     'tags': [CompilerBuild(Compiler.all, Build_type.mpiandsmp)]
    }
]

# These always have to be skipped else the suite hangs 
hanging_tests = [   
    # TODO(Alex) Issue #36 chargedensityplot hangs when running with np > 1 cores
    {'name':'groundstate-GGA_PBE-electronic_structure-Al',   
     'comment':'chargedensityplot hangs at an allgatherv call',
     'tags': [CompilerBuild(Compiler.all, Build_type.mpiandsmp)]}
]


def set_skipped_tests(executable:str, incl_failing_tests:bool) -> list:
    """
    Generate full list of tests to skip from the test suite

    :param incl_failing_tests: bool. If true, include failing tests in the suite,
                               hence exclude from the list of skipped tests

    :return tests_to_skip:     list of tests to skip, where each entry is a dict 
    """

    compiler = get_compiler_type()
    build_type = build_type_enum_map[executable]
    tests_to_skip = []

    if not incl_failing_tests:
        for test in failing_tests:
            for tag in test['tags']:
                if (tag.compiler == compiler or tag.compiler == Compiler.all) \
                and (tag.build == build_type or tag.build == Build_type.all):
                    tests_to_skip.append(test)
    
    # Always include hanging tests in tests to skip, else the test suite 
    # will hang 
    for test in hanging_tests:
        for tag in test['tags']:
            if (tag.compiler == compiler or tag.compiler == Compiler.all) \
            and (tag.build == build_type or tag.build == Build_type.all):
                tests_to_skip.append(test)

    return tests_to_skip
