# mkl.cmake
# This CMake script is responsible for finding and configuring Intel's Math Kernel Library (MKL).
# It allows users to enable MKL support in their project by setting the MKL option to ON.
# When MKL is found, the script will store the library and include directories for further use in the project.

# Option to enable the use of Intel MKL.
# This allows the user to decide whether to link against MKL in their build.
option(MKL "Use Intel MKL" OFF)

# Option for SCALAPACK.
# Check if the option for SCALAPACK (SCALAPACK) is already defined.
# If not, define the option to allow users to enable or disable SCALAPACK support.
if (NOT DEFINED SCALAPACK)
    option(SCALAPACK "Compile with SCALAPACK support" OFF)
endif()

# Check if the MKL option is enabled.
# Proceed with MKL-related configuration only if the user has set MKL to ON.
if(MKL)
    # Set ENABLE_SCALAPACK based on the value of SCALAPACK.
    # This will determine if SCALAPACK should be enabled when using MKL.
    set(ENABLE_SCALAPACK ${SCALAPACK})

    # Display the value of SCALAPACK in the configuration output.
    message(STATUS "SCALAPACK support with MKL: ${SCALAPACK}")

    # Here in case we want serial version we need to search for sequential
    if (NOT ${OMP})
        set(MKL_THREADING "sequential")
    endif()

    #Set the approapiate interface
    set(MKL_INTERFACE "lp64")

    # Attempt to find the MKL package using CMake's find_package command.
    # The CONFIG keyword ensures that CMake looks for a pre-configured MKL package.
    find_package(MKL CONFIG REQUIRED)

    # If MKL is found, log the libraries and the directory path.
    # These messages provide useful information about the MKL configuration during the CMake process.
    message(STATUS "MKL found (DIR): ${MKL_DIR}")
endif()
