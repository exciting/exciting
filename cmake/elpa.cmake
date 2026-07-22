# elpa.cmake
#
# This CMake script enables or disables ELPA support based on the user-defined
# option ELPA. If ELPA support is enabled, it finds the required ELPA 
# libraries and include files and sets the necessary compilation and linking parameters.

# Option to enable or disable ELPA support
option(ELPA "Enables ELPA support" OFF)

if (ELPA AND NOT SCALAPACK)
    message(FATAL_ERROR "ELPA requires SCALAPACK")
endif()

# Check if ELPA support is enabled
if(ELPA)

    # If enabled, display a status message
    message(STATUS "ELPA support : ON")
    
    # Allow user to override from command line
    set(ELPA_ROOT "" CACHE PATH "Path to the ELPA installation prefix")
    set(ELPA_OMP FALSE)

    # Trying to find elpa library
    if (OMP)
        find_library(ELPA_LIB NAMES elpa_openmp HINTS ${ELPA_ROOT}/lib/ ${ELPA_ROOT}/lib64/)
        if (NOT ELPA_LIB)
            message(STATUS "ELPA with OpenMP support not found. Trying ELPA without OpenMP")
            find_library(ELPA_LIB NAMES elpa HINTS ${ELPA_ROOT}/lib/ ${ELPA_ROOT}/lib64/)
        else()
            set(ELPA_OMP TRUE)
        endif()
    else()
        find_library(ELPA_LIB NAMES elpa HINTS ${ELPA_ROOT}/lib/ ${ELPA_ROOT}/lib64/)
    endif()

    if(NOT ELPA_LIB)
        message(FATAL_ERROR "ELPA not found. Please provide ELPA_ROOT to the ELPA install directory.")
    endif()

    # Get directory
    get_filename_component(ELPA_LIBDIR ${ELPA_LIB} DIRECTORY)

    # If ELPA is found, display the paths
    message(STATUS "ELPA library found at: ${ELPA_LIBDIR}")

    if (ELPA_OMP)
        find_file(ELPA_INFO_FILE NAMES elpa_openmp.pc HINTS ${ELPA_LIBDIR}/pkgconfig/)
    else()
        find_file(ELPA_INFO_FILE NAMES elpa.pc HINTS ${ELPA_LIBDIR}/pkgconfig/)
    endif()
    
    if(ELPA_INFO_FILE)
        # Read the content of elpa.pc file
        file(READ ${ELPA_INFO_FILE} ELPA_INFO_CONTENTS)

        # Get the prefix
        string(REGEX MATCH "prefix=([^ \n]+)" MATCH_PREFIX ${ELPA_INFO_CONTENTS})
        set(ELPA_PREFIX_DIR ${CMAKE_MATCH_1})

        # Parse the version from the file
        string(REGEX MATCH "Version: ([0-9]+)\\.([0-9]+)\\.([0-9]+)" VERSION_MATCH ${ELPA_INFO_CONTENTS})
        set(ELPA_MAJOR_VERSION ${CMAKE_MATCH_1})
        set(ELPA_MINOR_VERSION ${CMAKE_MATCH_2})
        set(ELPA_PATCH_VERSION ${CMAKE_MATCH_3})
        set(ELPA_VERSION "${ELPA_MAJOR_VERSION}.${ELPA_MINOR_VERSION}.${ELPA_PATCH_VERSION}")
        message(STATUS "ELPA version: ${ELPA_VERSION}")
    else()
        # If no pc file, try to detect version from directory name
        set(ELPA_VERSION "unknown")
        message(WARNING "Could not find elpa.pc file, version detection may be inaccurate")
    endif()

    # Find header file path - updated for your directory structure
    if(ELPA_OMP)
        set(ELPA_INCLUDE_PATH "${ELPA_PREFIX_DIR}/include/elpa_openmp-${ELPA_VERSION}")
        set(ELPA_MODULES_PATH "${ELPA_PREFIX_DIR}/include/elpa_openmp-${ELPA_VERSION}/modules")
    else()
        set(ELPA_INCLUDE_PATH "${ELPA_PREFIX_DIR}/include/elpa-${ELPA_VERSION}")
        set(ELPA_MODULES_PATH "${ELPA_PREFIX_DIR}/include/elpa-${ELPA_VERSION}/modules")
    endif()

    # Check if paths exist
    if(NOT EXISTS "${ELPA_INCLUDE_PATH}")
        message(FATAL_ERROR "ELPA include path not found: ${ELPA_INCLUDE_PATH}")
    endif()
    if(NOT EXISTS "${ELPA_MODULES_PATH}")
        message(FATAL_ERROR "ELPA modules path not found: ${ELPA_MODULES_PATH}")
    endif()

    # Include directories
    include_directories(${ELPA_INCLUDE_PATH})
    include_directories(${ELPA_INCLUDE_PATH}/elpa)
    include_directories(${ELPA_MODULES_PATH})

    message(STATUS "ELPA include directories: ${ELPA_INCLUDE_PATH};${ELPA_MODULES_PATH}")

    # Set libraries
    link_directories(${ELPA_LIBDIR})
    if(ELPA_OMP)
        set(ELPA_LIBS "-lelpa_openmp")
    else()
        set(ELPA_LIBS "-lelpa")
    endif()

    message(STATUS "ELPA libraries: ${ELPA_LIBS}")

    # Add compiler preprocessor flag for conditional compilation
    add_compile_definitions(_ELPA_)
    
else()
    # If ELPA support is not enabled, display a status message
    message(STATUS "ELPA support : OFF")
endif()
