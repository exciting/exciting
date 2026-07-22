# mpi.cmake
# This file handles the inclusion and configuration of MPI (Message Passing Interface) support
# in a CMake-based project. If the MPI option is enabled, it attempts to find and configure
# MPI compilers for C++, C, and Fortran. It sets the necessary flags and includes directories,
# as well as a definition to indicate MPI support in the code.
# Note that one can add small program test following the basic example here to
# test for specfic MPI characteristics support

# Option to enable MPI build
option(MPI "Enable MPI support" ON)

# Option to enable MPI checks.
# Disable if the system does not allow to launch MPI programs
# on login nodes.
option(RUN_MPI_CHECKS "Compile and run a minimal MPI program.
    Set to OFF if the system does not allow to launch MPI programs
    on login nodes." ON)

if(MPI)

    # Find the MPI package (required for compiling with MPI support).
    find_package(MPI REQUIRED)
    message(STATUS "MPI support: ON")

    if(NOT MPI_Fortran_FOUND)
        message(FATAL_ERROR "No MPI Fortran compiler/library was found.")
    endif()

    if(NOT MPI_Fortran_HAVE_F08_MODULE)
        message(FATAL_ERROR "mpi_f08 is not available in this MPI implementation.")
    endif()

    if(RUN_MPI_CHECKS)
        message(STATUS "Testing the MPI configuration")

        set(_mpi_dir ${CMAKE_BINARY_DIR}/config_tests)
        file(MAKE_DIRECTORY ${_mpi_dir})
        set(MPI_CONFIG_TEST_LOG "${_mpi_dir}/mpi_test_log.txt")
        file(WRITE "${MPI_CONFIG_TEST_LOG}" "MPI configuration check\n\n")

        # 1. Build the real source file (cmake/checks/test_mpi_check.f90)
        try_compile(MPI_TEST_BUILDS
            ${_mpi_dir}
            SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/checks/test_mpi_check.f90
            LINK_LIBRARIES MPI::MPI_Fortran
            COPY_FILE ${_mpi_dir}/test_mpi_check
            OUTPUT_VARIABLE _mpi_build_log
        )
        file(APPEND "${MPI_CONFIG_TEST_LOG}" "== build ==\n${_mpi_build_log}\n")

        if(NOT MPI_TEST_BUILDS)
            message(FATAL_ERROR "MPI test failed to compile. See ${MPI_CONFIG_TEST_LOG}")
        endif()

        # 2. Run the test. Execute the binary directly first: a single-rank
        #    program works as a singleton and this is required on machines where
        #    the MPI launcher may not run on the login/build node (e.g. Cray,
        #    where MPIEXEC_EXECUTABLE is srun). Only if the direct run fails
        #    (e.g. some Spack setups that mandate a launcher) fall back to the
        #    launcher FindMPI discovered (mpiexec/mpirun/srun/...).
        execute_process(
            COMMAND ${_mpi_dir}/test_mpi_check
            WORKING_DIRECTORY ${_mpi_dir}
            RESULT_VARIABLE _mpi_run_result
            OUTPUT_VARIABLE _mpi_run_log
            ERROR_VARIABLE  _mpi_run_log
        )
        file(APPEND "${MPI_CONFIG_TEST_LOG}" "== run (direct) ==\n${_mpi_run_log}\n")

        # If using spack, depending on the configuration one might require mpirun
        # for the code to run
        if(NOT _mpi_run_result EQUAL 0)
            execute_process(
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                        ${MPIEXEC_PREFLAGS}
                        ${_mpi_dir}/test_mpi_check
                        ${MPIEXEC_POSTFLAGS}
                WORKING_DIRECTORY ${_mpi_dir}
                RESULT_VARIABLE _mpi_run_result
                OUTPUT_VARIABLE _mpi_run_log
                ERROR_VARIABLE  _mpi_run_log
            )
            file(APPEND "${MPI_CONFIG_TEST_LOG}" "== run (${MPIEXEC_EXECUTABLE}) ==\n${_mpi_run_log}\n")
        endif()

        if(NOT _mpi_run_result EQUAL 0)
            message(FATAL_ERROR "MPI test failed to run. See ${MPI_CONFIG_TEST_LOG}")
        endif()

        message(STATUS "MPI test passed during configuration.")

        # =========================================================================
        # Check for MPI 4.0 Large Count Support & Wrapper Stability
        # =========================================================================
        message(STATUS "Testing MPI 4.0 large count support ...")

        # 1. Build the external tests that check IN_PLACE and 64-bit Truncation
        try_compile(MPI4_IN_PLACE_BUILDS
            ${_mpi_dir}
            SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/checks/test_mpi4_in_place.f90
            LINK_LIBRARIES MPI::MPI_Fortran
            COPY_FILE ${_mpi_dir}/test_mpi4_in_place
            OUTPUT_VARIABLE _mpi4_in_place_build_log
        )
        file(APPEND "${MPI_CONFIG_TEST_LOG}" "MPI-4 IN_PLACE Compile Log:\n${_mpi4_in_place_build_log}\n")

        try_compile(MPI4_TRUNC_BUILDS
            ${_mpi_dir}
            SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/checks/test_mpi4_truncation.f90
            LINK_LIBRARIES MPI::MPI_Fortran
            COPY_FILE ${_mpi_dir}/test_mpi4_truncation
            OUTPUT_VARIABLE _mpi4_trunc_build_log
        )
        file(APPEND "${MPI_CONFIG_TEST_LOG}" "MPI-4 Truncation Compile Log:\n${_mpi4_trunc_build_log}\n")

        # 2. Evaluate the Results
        if(MPI4_IN_PLACE_BUILDS AND MPI4_TRUNC_BUILDS)

            # Check 1: IN_PLACE wrapper bug
            execute_process(
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${_mpi_dir}/test_mpi4_in_place
                WORKING_DIRECTORY ${_mpi_dir}
                RESULT_VARIABLE MPI_IN_PLACE_RESULT
                OUTPUT_QUIET ERROR_QUIET TIMEOUT 15
            )

            # Check 2: 64-bit to 32-bit Truncation bug
            execute_process(
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${_mpi_dir}/test_mpi4_truncation
                WORKING_DIRECTORY ${_mpi_dir}
                RESULT_VARIABLE MPI_TRUNC_RESULT
                OUTPUT_QUIET ERROR_QUIET TIMEOUT 15
            )

            if(MPI_IN_PLACE_RESULT EQUAL 0 AND MPI_TRUNC_RESULT EQUAL 0)
                message(STATUS "MPI 4.0 Status: FULLY SUPPORTED AND STABLE")
                message(STATUS "Enabling native 64-bit counts (USE_MPI4_LARGE_COUNTS).")
                add_compile_definitions(USE_MPI4_LARGE_COUNTS)
            else()
                message(WARNING "MPI 4.0 Status: SUPPORTED BUT UNSTABLE")
                message(WARNING "The MPI library compiled successfully but failed runtime safety checks:")
                message(WARNING "  IN_PLACE Exit Code: ${MPI_IN_PLACE_RESULT}")
                message(WARNING "  TRUNCATION Exit Code: ${MPI_TRUNC_RESULT}")
                message(STATUS "Falling back to the safe 32-bit interface. USE_MPI4_LARGE_COUNTS will NOT be enabled.")
                file(APPEND "${MPI_CONFIG_TEST_LOG}" "MPI-4 Check Failed. IN_PLACE: ${MPI_IN_PLACE_RESULT}, Trunc: ${MPI_TRUNC_RESULT}\n")
            endif()

        else()
            message(STATUS "MPI 4.0 Status: NOT SUPPORTED")
            message(STATUS "The MPI compiler lacks large count interfaces. Falling back to the safe 32-bit interface.")
        endif()

    else()
        message(STATUS "MPI tests skipped (RUN_MPI_CHECKS=OFF).")
        message(STATUS "USE_MPI4_LARGE_COUNTS will NOT be enabled.")
    endif()

    # Apply global MPI settings only after the checks succeed (if ran)
    include_directories(SYSTEM ${MPI_Fortran_INCLUDE_DIRS})
    add_compile_definitions(MPI)

    # Libraries linked into the exciting targets (consumed in src/CMakeLists.txt).
    # Link the imported target so that all MPI Fortran libraries and flags are
    # pulled in: plain compilers such as gfortran do not auto-link MPI (unlike
    # the Intel/Cray compiler wrappers), which otherwise causes undefined
    # references to mpi_* symbols at link time. The explicit library lists are
    # kept as well so that scalapack.cmake can still detect the MPI flavour
    # (openmpi/mpich) from the library paths.
    set(MPI_LIBS MPI::MPI_Fortran ${MPI_Fortran_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})

    # Get the vendor
    if(NOT MPIEXEC_EXECUTABLE)
        message(WARNING "MPIEXEC_EXECUTABLE is not defined, so the MPI vendor and version cannot be "
                         "detected. Known-buggy MPI versions will therefore not be caught automatically.")
    else()
        execute_process(
            COMMAND ${MPIEXEC_EXECUTABLE} --version
            OUTPUT_VARIABLE MPI_VERSION_OUTPUT
            ERROR_VARIABLE MPI_VERSION_OUTPUT # Some vendors output version info to stderr
            RESULT_VARIABLE MPI_VERSION_RESULT
        )

        if(NOT MPI_VERSION_RESULT EQUAL 0)
            message(WARNING "Failed to run ${MPIEXEC_EXECUTABLE} --version")
        endif()

        string(TOLOWER "${MPI_VERSION_OUTPUT}" MPI_VERSION_OUTPUT_LOWER)

        if(MPI_VERSION_OUTPUT_LOWER MATCHES "intel")
            set(MPI_VENDOR "Intel")
            # Check for buggy versions of Intel MPI.
            # Different Intel MPI releases use different banner formats, e.g.:
            #   "Intel(R) MPI Library for Linux* OS, Version 2019 Update 6 Build ..."
            #   "Intel(R) MPI Library 2021.18 for Linux* OS"
            # The latter has no "Version" keyword, so instead of anchoring on
            # that word we look for a bare year-based version number
            # (####.## or ####.##.##) anywhere in the banner.
            if(MPI_VERSION_OUTPUT_LOWER MATCHES "([0-9][0-9][0-9][0-9]\\.[0-9]+(\\.[0-9]+)?)")
                set(INTEL_MPI_VERSION "${CMAKE_MATCH_1}")
                message(STATUS "Intel MPI Version detected: ${INTEL_MPI_VERSION}")

                # The known bug affects the whole 2021.18.x release line, not
                # just 2021.18.0, so match on the major.minor prefix and
                # ignore the patch component.
                if(INTEL_MPI_VERSION MATCHES "^2021\\.18(\\.|$)")
                    message(FATAL_ERROR "Intel MPI version ${INTEL_MPI_VERSION} is known to be buggy and is unsupported. Please upgrade or downgrade your MPI toolkit.")
                endif()
            else()
                message(WARNING "Could not parse Intel MPI version number.")
            endif()
        elseif(MPI_VERSION_OUTPUT_LOWER MATCHES "open-mpi" OR MPI_VERSION_OUTPUT_LOWER MATCHES "open mpi" OR MPI_VERSION_OUTPUT_LOWER MATCHES "openmpi")
            set(MPI_VENDOR "OpenMPI")
        # MVAPICH and Cray's MPICH-based wrappers can also report "mpich" in
        # their --version banner, so they must be checked before the generic
        # MPICH match below, or they would be misdetected as plain MPICH.
        elseif(MPI_VERSION_OUTPUT_LOWER MATCHES "mvapich")
            set(MPI_VENDOR "MVAPICH")
        elseif(MPI_VERSION_OUTPUT_LOWER MATCHES "cray")
            set(MPI_VENDOR "Cray")
        elseif(MPI_VERSION_OUTPUT_LOWER MATCHES "mpich")
            set(MPI_VENDOR "MPICH")

            if(MPI_VERSION_OUTPUT_LOWER MATCHES "version:[ \t]+([0-9]+(\\.[0-9]+)+)")
                set(MPI_VERSION_STRING "${CMAKE_MATCH_1}")

                # Convert to comparable number (major.minor)
                string(REGEX MATCH "^([0-9]+)\\.([0-9]+)" _ "${MPI_VERSION_STRING}")
                set(MPI_VERSION_MAJOR "${CMAKE_MATCH_1}")
                set(MPI_VERSION_MINOR "${CMAKE_MATCH_2}")

                # Convert to numeric value: major.minor -> major*100 + minor (safe comparison)
                math(EXPR MPI_VERSION_NUM "${MPI_VERSION_MAJOR} * 100 + ${MPI_VERSION_MINOR}")

                # MPI 5.0 -> 500
                if(MPI_VERSION_NUM GREATER_EQUAL 500)
                    add_compile_definitions(MPI_F_SYNC_REG_DEPRECATED)
                    message(STATUS "MPI_F_SYNC_REG is depreacted")
                endif()
            endif()

        else()
            set(MPI_VENDOR "Unknown")
        endif()

        message(STATUS "Detected MPI Vendor: ${MPI_VENDOR}")

    endif()

else()
    # Notify that MPI support is disabled.
    message(STATUS "MPI support: OFF")
endif()
