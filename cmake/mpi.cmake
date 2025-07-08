# mpi.cmake
# This file handles the inclusion and configuration of MPI (Message Passing Interface) support
# in a CMake-based project. If the MPI option is enabled, it attempts to find and configure
# MPI compilers for C++, C, and Fortran. It sets the necessary flags and includes directories,
# as well as a definition to indicate MPI support in the code.
# Note that one can add small program test following the basic example here to 
# test for specfic MPI characteristics support

# Option controling in we are going to do an MPI compilation
option(MPI "Enable MPI support" ON)

if(MPI)
    # Find the MPI package, which is required for compiling with MPI support.
    find_package(MPI REQUIRED)
    message(STATUS "MPI support: ON")

    # Check if MPI C++ compiler is found, if not, terminate with a fatal error.
    if(NOT MPI_CXX_FOUND)
        message(FATAL_ERROR "No MPI compilers were found (C++).")
    endif()

    # Include the MPI C++ header files.
    include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})

    # Check if MPI C compiler is found, if not, terminate with a fatal error.
    if(NOT MPI_C_FOUND)
        message(FATAL_ERROR "No MPI compilers were found (C).")
    endif()

    # Check if MPI Fortran compiler is found, if not, terminate with a fatal error.
    if(NOT MPI_Fortran_FOUND)
        message(FATAL_ERROR "No MPI compilers were found (Fortran).")
    endif()

    # For fortran ensure the Fortran 2008 modules to be available
    if(NOT MPI_Fortran_HAVE_F08_MODULE)
    	message(FATAL_ERROR "mpi_f08 is not available in this MPI implementation.")
    endif()

    # Include the MPI Fortran header files.
    include_directories(SYSTEM ${MPI_Fortran_INCLUDE_PATH})
    add_compile_definitions(MPI)

    message(STATUS "Testing the MPI configuration")

    # Set MPI libraries for linking.
    set(MPI_LIBS "${MPI_CXX_LIBRARIES};${MPI_C_LIBRARIES};${MPI_Fortran_LIBRARIES}")

    # Create a simple Fortran MPI test program.
    file(WRITE ${CMAKE_BINARY_DIR}/config_tests/test_mpi_check.f90 "
    program test_mpi_check
        use mpi
        implicit none

        integer :: rank, size, ierr

        call MPI_Init(ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

        if (size < 1) then
            print *, 'MPI setup error: size is less than 1.'
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        else
            print *, 'MPI test passed: Hello from process', rank, 'of', size
        end if

        call MPI_Finalize(ierr)
    end program test_mpi_check
    ")

    # Create a log file for the config test 
    set(MPI_CONFIG_TEST_LOG "${CMAKE_BINARY_DIR}/config_tests/mpi_test_log.txt")
    file(WRITE "${MPI_CONFIG_TEST_LOG}" "Compilation Log:\n")
    set(MPI_CONFIG_TEST_LOG_STR "")

    # Try to compile and run the Fortran test program during configuration.
    execute_process(
            COMMAND ${MPI_Fortran_COMPILER} -o test_mpi_check test_mpi_check.f90
	    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/config_tests
	    RESULT_VARIABLE MPI_TEST_COMPILE_RESULT
	    OUTPUT_VARIABLE MPI_CONFIG_TEST_LOG_STR
	    ERROR_VARIABLE  MPI_CONFIG_TEST_LOG_STR
    )

    # Save log
    file(APPEND "${MPI_CONFIG_TEST_LOG}" "${MPI_CONFIG_TEST_LOG_STR}\n")

    execute_process(
            COMMAND ./test_mpi_check
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/config_tests
	    RESULT_VARIABLE MPI_TEST_RUN_RESULT
            OUTPUT_VARIABLE MPI_CONFIG_TEST_LOG_STR
            ERROR_VARIABLE  MPI_CONFIG_TEST_LOG_STR
    )

    # If using spack, depending on the configuration one might require mpirun
    # for the code to run
    if (MPI_TEST_COMPILE_RESULT EQUAL 0 AND NOT MPI_TEST_RUN_RESULT EQUAL 0)
        execute_process(
                COMMAND mpirun -np 1 ./test_mpi_check
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/config_tests
                RESULT_VARIABLE MPI_TEST_RUN_RESULT
                OUTPUT_VARIABLE MPI_CONFIG_TEST_LOG_STR
                ERROR_VARIABLE  MPI_CONFIG_TEST_LOG_STR
        )
    endif()

    # Save log
    file(APPEND "${MPI_CONFIG_TEST_LOG}" "${MPI_CONFIG_TEST_LOG_STR}\n")

    # Check if the test compiled and ran successfully.
    if(NOT MPI_TEST_COMPILE_RESULT EQUAL 0)
	message(FATAL_ERROR "MPI test failed to compile. Please check your MPI setup.")
    elseif(NOT MPI_TEST_RUN_RESULT EQUAL 0)
	message(FATAL_ERROR "MPI test failed to run. Please check your MPI setup.")
    else()
        message(STATUS "MPI test passed during configuration.")
    endif()

else()
    # Notify that MPI support is disabled.
    message(STATUS "MPI support: OFF")
endif()

