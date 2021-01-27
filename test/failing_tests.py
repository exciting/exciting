"""
 Failing and problem tests.

 Each skipped test should have a TODO(Name) assigned to it here, which
 is also logged on the gitlab issue tracker.

 Each list entry should be a dictionary.
"""

# Partially failing tests

# TODO(Ronaldo) Issue #34. groundstate-GGA_PBE-tddft-LiF.
#   MPI calculations result in numerous files containing numbers that differ from serial calculations

# TODO(Alex) Issue #35. groundstate-LDA_PW-gw-Si.
#   MPI GW calculations with np > 1 does not produce EPS00_GW.OUT 

skipped_tests_mpismp = [
    # TODO(Alex) Issue #36 chargedensityplot hangs when running with np > 1 cores
    {'name':'groundstate-GGA_PBE-electronic_structure-Al',   
     'comment':'Hanging at an allgatherv call'}
]

skipped_tests_serial = []

skipped_tests_mpi = []

# Only edit if a new build type is added to the test suite 
skipped_tests = {'excitingser':skipped_tests_serial,
                 'excitingmpi':skipped_tests_mpi,
                 'excitingmpismp':skipped_tests_mpismp
                }
    
