# scalapack.cmake
# This CMake script is designed to locate and configure the SCALAPACK library
# and its interface with Intel MKL if available. The script will check for
# the presence of MKL and set the appropriate include directories and library
# paths for SCALAPACK based on the user's configuration.

# Set the root directory for the SCALAPACK installation, allowing for non-standard locations.
# If the user has a custom location for SCALAPACK, they can set SCALAPACK_ROOT.
set(SCALAPACK_ROOT "None" CACHE STRING "Root directory for non-standard locations for the scalapack installation")

# Option to enable or disable SCALAPACK support in the build.
option(SCALAPACK "Compile with SCALAPACK support" OFF)
option(CRAY_SCALAPACK "Tells the compiler we are using SCALAPACK from the Cray libsci" OFF)
message(STATUS "ScaLAPACK support: ${SCALAPACK}")

# Check if MKL (Math Kernel Library) is being used for SCALAPACK.
if(MKL AND SCALAPACK AND NOT CRAY_SCALAPACK)
    # If MKL is specified but not found, search for it.
    if (NOT MKL_FOUND)
	set(MKL_INTERFACE "lp64")
	if (NOT ${OMP})
            set(MKL_THREADING "sequential")
        endif()
	find_package(MKL CONFIG REQUIRED)
        message(STATUS "MKL found (DIR): ${MKL_DIR}")
    endif()

    # Indicate that SCALAPACK is being used with MKL.
    set(SCALAPACK_FOUND True)

    # Include the necessary directories from the MKL installation.
    include_directories(${MKL_ROOT}/include/)

    # Find the SCALAPACK library provided by MKL.
    find_library(SCALAPACK_LIB NAMES mkl_scalapack_lp64 HINTS ${MKL_ROOT}/lib)
    message(STATUS "Using SCALAPACK from Intel MKL: ${SCALAPACK_LIB}")

    # Define a preprocessor macro to signal that SCALAPACK is enabled.
    add_compile_definitions(SCAL)
endif()

# If MKL is not being used but SCALAPACK support is enabled:
if(NOT MKL AND SCALAPACK AND NOT CRAY_SCALAPACK)
    # Attempt to find the SCALAPACK library manually in the specified directory.
    find_library(SCALAPACK_LIB NAMES scalapack scalapack-mpi HINTS ${SCALAPACK_ROOT}/lib/)

    # If not found, check for more specific names based on the MPI flavor.
    if (NOT SCALAPACK_LIB)
        # Check for different MPI flavors and attempt to find the appropriate SCALAPACK version.
        # If the MPI installation contains "openmpi" in the path, look for scalapack-openmpi.
        # If it contains "mpich", look for scalapack-mpich or scalapack-mpich2.
        string(TOLOWER "${MPI_LIBS}" MPI_LIBS_LOWER)
        if(MPI_LIBS_LOWER MATCHES "openmpi")
            find_library(SCALAPACK_LIB NAMES scalapack-openmpi HINTS ${SCALAPACK_ROOT}/lib/)
        elseif(MPI_LIBS_LOWER MATCHES "mpich")
            find_library(SCALAPACK_LIB NAMES scalapack-mpich scalapack-mpich2 HINTS ${SCALAPACK_ROOT}/lib/)
        endif()
    endif()

    # If the SCALAPACK library is found, print a success message.
    if (SCALAPACK_LIB)
        message(STATUS "Found SCALAPACK: ${SCALAPACK_LIB}")
        set(SCALAPACK_FOUND True)
    else()
        # If not found, display an error message and stop the configuration.
        message(FATAL_ERROR "SCALAPACK not found, adjust SCALAPACK_ROOT")
    endif()

    # Define a preprocessor macro to signal that SCALAPACK is enabled.
    add_compile_definitions(SCAL)

endif()

if(NOT MKL AND SCALAPACK AND CRAY_SCALAPACK)
    # Find the current BLAS libraries
    string(TOLOWER "${BLAS_LIBRARIES}" BLAS_LIBS_LOWER)
    if (NOT BLAS_LIBS_LOWER MATCHES "libsci")
        message(FATAL_ERROR "For using libsci SCALAPACK, the linear algebra needs to be done through libsci")
    endif()
    message(STATUS "Found SCALAPACK: ${BLAS_LIBRARIES}")
    set(SCALAPACK_FOUND True)
    add_compile_definitions(SCAL)
endif()

