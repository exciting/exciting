# mkl.cmake
# This CMake script is responsible for finding and configuring Intel's Math Kernel Library (MKL).
# It allows users to enable MKL support in their project by setting the MKL option to ON.
# When MKL is found, the script will store the library and include directories for further use in the project.

# Option to enable the use of Intel MKL.
# This allows the user to decide whether to link against MKL in their build.
option(MKL "Use Intel MKL" OFF)

# Check if the MKL option is enabled.
# Proceed with MKL-related configuration only if the user has set MKL to ON.
if(MKL)
    # Here in case we want serial version we need to search for sequential
    if (NOT OMP)
        set(MKL_THREADING "sequential")
    endif()

    #Set the approapiate interface
    set(MKL_INTERFACE "lp64")

    # Attempt to find the MKL package using CMake's find_package command.
    # The CONFIG keyword ensures that CMake looks for a pre-configured MKL package.
    find_package(MKL CONFIG)

    # If MKL is not found, try looking in the default MKLROOT path
    # This solve issues in old MKL installations that do not 
    # load the proper path to CMAKE_PREFIX_PATH.
    if (NOT MKL_FOUND AND DEFINED ENV{MKLROOT})
        set(MKL_ROOT_DIR "$ENV{MKLROOT}/lib/cmake/mkl")
        if (EXISTS "${MKL_ROOT_DIR}")
            list(APPEND CMAKE_PREFIX_PATH "${MKL_ROOT_DIR}")
            find_package(MKL CONFIG REQUIRED)
        else()
	    message(FATAL_ERROR "CMake configuration file for MKL cannot be found. Please check your MKL configuration.")
        endif()
    endif()

    # Print information regarding MKL 
    if (MKL_FOUND)
        message(STATUS "MKL found (DIR): ${MKL_DIR}")
    else()
        message(FATAL_ERROR "Intel MKL could not be found. Please check your MKL installation.")
    endif()


endif()
