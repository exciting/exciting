# fftw3.cmake
# This CMake script is designed to locate and configure the FFTW3 library
# and its interface with Intel MKL if available. The script will check for
# the presence of MKL and set the appropriate include directories and library 
# paths for FFTW3 and FFTW3F based on the user's configuration.

# Set the root directory for the FFTW3 installation, allowing for non-standard locations.
set(FFTW3_ROOT "None" CACHE STRING "Root directory for non-standard locations for the FFTW3 installation")

# Check if MKL (Math Kernel Library) is being used
if(MKL)
    # If MKL is specified but not found, search for it
    if (NOT MKL_FOUND)
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

    # Indicate that FFTW3 is being used with MKL
    set(FFTW3_FOUND True)
    
    # Include the necessary directories from the MKL installation
    include_directories(${MKL_ROOT}/include/fftw/)
    include_directories(${MKL_ROOT}/include/)
    message(STATUS "Using FFTW3 interface to MKL FFT")
    
    # Define preprocessor macro for FFTW3 interface
    add_compile_definitions(FFTW3_INTERFACE FFTW)

else()   
    # If MKL is not used, search for FFTW3 libraries and headers
    find_library(FFTW3_LIBRARIES NAMES fftw3 HINTS ${FFTW3_ROOT}/lib)
    find_library(FFTW3F_LIBRARIES NAMES fftw3f HINTS ${FFTW3_ROOT}/lib)

    # If OpenMP is required, search for the corresponding OpenMP libraries
    if (OMP)
        find_library(FFTW3_LIBRARIES_OMP NAMES fftw3_omp HINTS ${FFTW3_ROOT}/lib)
        find_library(FFTW3F_LIBRARIES_OMP NAMES fftw3f_omp HINTS ${FFTW3_ROOT}/lib)
    endif()

    # Locate the FFTW3 header file
    find_path(FFTW3_INCLUDE_DIR NAMES "fftw3.h" HINTS ${FFTW3_ROOT}/include)
    
    # Check if the FFTW3 libraries and headers were found
    if (NOT FFTW3_LIBRARIES OR NOT FFTW3_INCLUDE_DIR)
        message(FATAL_ERROR "FFTW not found")
    else()
        set(FFTW3_FOUND True)
    endif()

    # Log the inclusion of FFTW3 directories
    message(STATUS "FFTW3 found. Include directories: ${FFTW3_INCLUDE_DIR}")
    list(APPEND FFTW3_LIBRARIES ${FFTW3F_LIBRARIES})
    
    # Append OpenMP libraries if required
    if (OMP)
        list(APPEND FFTW3_LIBRARIES ${FFTW3_LIBRARIES_OMP})
        list(APPEND FFTW3_LIBRARIES ${FFTW3F_LIBRARIES_OMP})
    endif()

    # Log the libraries found for FFTW3
    message(STATUS "FFTW3 found. Libraries: ${FFTW3_LIBRARIES}")

    # Include the FFTW3 header directory
    include_directories(${FFTW3_INCLUDE_DIR})

    # Define preprocessor macro for FFTW3 interface
    add_compile_definitions(FFTW3_INTERFACE)
endif()
