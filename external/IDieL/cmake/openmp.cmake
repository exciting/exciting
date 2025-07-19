# openmp.cmake
#
# This CMake script configures OpenMP support for a Fortran project.
# It provides options to enable or disable OpenMP and creates a test
# source file to verify the availability of OpenMP offload capabilities
# on the target architecture. If OpenMP is enabled, the necessary
# compile definitions and flags are set, and a simple test program
# is generated and compiled to check for compatibility.
#
# Usage:
# - Include this file in your main CMakeLists.txt to enable OpenMP support.
# - The script will generate a Fortran file and attempt to compile it
#   during the configuration phase.

# Option to enable or disable OpenMP support
option(OMP "Enable OpenMP support" ON)

# Check if the OpenMP support option is enabled
if(OMP)

    # Display a status message indicating that OpenMP support is ON
    message(STATUS "OpenMP support : ON")

    # Find the OpenMP package; this will fail if OpenMP is not available
    find_package(OpenMP REQUIRED)

    # If enabled, define the USEOMP compilation flag
    add_compile_definitions(USEOMP)

    # Check the compiler
    if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        # Check if it is not IntelLLVM
        if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
	   message(STATUS "Using Classic Intel Fortran compiler. Overwriting OpenMP Fortran flags by the appropiate ones, i.e. not ifx.")
           # Overwrite OpenMP Fortran flags as needed
           set(OpenMP_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -qopenmp -diag-disable=10448")
        else()
           message(STATUS "Using Intel LLVM Fortran compiler. OpenMP flags are not modified.")
        endif()
    endif()

    # Create a Fortran source file for testing OpenMP offload capabilities
    file(WRITE ${CMAKE_BINARY_DIR}/config_tests/test_omp_offload.f90 "
    program main

        use iso_fortran_env,         only: i32 => int32, r64 => real64
        use iso_c_binding
        use omp_lib

        implicit none

        integer(i32), parameter   :: n = 4000_i32
        complex(r64), parameter   :: zzero = cmplx(0.0, 0.0, kind=r64)
        complex(r64), parameter   :: zone  = cmplx(1.0, 0.0, kind=r64)

        complex(r64), pointer, contiguous :: A(:,:)
        type(c_ptr) :: A_cptr

        complex(r64), allocatable :: B(:,:)
        complex(r64), allocatable :: C(:,:)
        integer(i32)              :: i, j
        real(r64)                :: start, finish

        call cpu_time(start)

        ! Imagine B is given by an external routine in the CPU
        allocate(B(n,n), source=cmplx(42.0, 42.0, kind=r64))
        !$omp target enter data map(to: B)

        allocate(C(n,n))
        !$omp target enter data map(alloc: C)

        A_cptr = omp_target_alloc(2 * c_sizeof(1.0_c_double) * n * n, omp_get_default_device())
        call c_f_pointer(A_cptr, A, int([n,n],kind=c_size_t))
        A(0:n-1,0:n-1) => A(:,:)

        !$omp target has_device_addr(A)
        !$omp teams distribute parallel do simd private(i, j) collapse(2)
        do i = 0, n-1
            do j = 0, n-1
                A(j,i) = fillA(i,j) + 1.0
                C(j,i) = A(j,i) + B(j,i)
            end do
        end do
        !$omp end teams distribute parallel do simd
        !$omp end target

        nullify(A)
        call omp_target_free(A_cptr, omp_get_default_device())
        !$omp target update from(C)
        !$omp target exit data map(delete: C)
        !$omp target exit data map(delete: B)

        call cpu_time(finish)

        write(*,*) A(1,1)
        print '(\"Time = \", f16.3, \" seconds.\")', finish - start

    contains

        pure complex(r64) function fillA(i,j)

            integer(i32), intent(in) :: i, j
            !$omp declare target(fillA)

            fillA = cmplx(0, 0, r64)

        end function fillA

    end program main")

    # Create a log file for the config test
    set(OMPOFFLOAD_CONFIG_TEST_LOG "${CMAKE_BINARY_DIR}/config_tests/ompoffload_test_log.txt")
    file(WRITE "${OMPOFFLOAD_CONFIG_TEST_LOG}" "Compilation Log:\n")
    set(OMPOFFLOAD_CONFIG_TEST_LOG_STR "")

    # Try to compile and run the Fortran test program during configuration.
    execute_process(
            COMMAND ${CMAKE_Fortran_COMPILER} ${OpenMP_Fortran_FLAGS} -o test_omp_offload test_omp_offload.f90
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/config_tests
            RESULT_VARIABLE OMPOFFLOAD_TEST_COMPILE_RESULT
            OUTPUT_VARIABLE OMPOFFLOAD_CONFIG_TEST_LOG_STR
            ERROR_VARIABLE  OMPOFFLOAD_CONFIG_TEST_LOG_STR
    )

    # Append the compilation log to the log file
    file(APPEND "${OMPOFFLOAD_CONFIG_TEST_LOG}" "${OMPOFFLOAD_CONFIG_TEST_LOG_STR}\n")

    # Check the compilation result and display a status message
    if (OMPOFFLOAD_TEST_COMPILE_RESULT)
        message(STATUS "OMP OFFLOAD: Supported")
    else()
        message(STATUS "OMP OFFLOAD: Unsupported")
    endif()

else()
    # Display a status message indicating that OpenMP support is OFF
    message(STATUS "OpenMP support : OFF")
endif()
