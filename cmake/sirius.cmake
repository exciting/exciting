# sirius.cmake
# SIRIUS is a domain-specific library for electronic structure calculations.
# It implements pseudopotential plane wave (PP-PW) and full potential
# linearized augmented plane wave (FP-LAPW) methods and is designed for
# GPU acceleration of popular community codes such as Exciting, Elk,
# and Quantum ESPRESSO. SIRIUS is written in C++17 with MPI, OpenMP,
# and CUDA/ROCm programming models.
# SIRIUS is organized as a collection of classes that abstract away the
# different building blocks of the DFT self-consistency cycle.

# Option to enable SIRIUS support in the project.
option(SIRIUS "Enables SIRIUS support in exciting" OFF)
# Display the current status of the SIRIUS support option.
message(STATUS "SIRIUS support : ${SIRIUS}")

# Check if SIRIUS support is enabled.
if (SIRIUS)

    # Set the root directory for the SIRIUS installation.
    # This allows specifying a non-standard installation location.
    set(SIRIUS_ROOT "None" CACHE STRING "Root directory for SIRIUS installation if not in the path")

    # Attempt to find the SIRIUS library using find_library.
    # It searches for the library named 'sirius' in the specified hints.
    find_library(SIRIUS_LIB NAMES sirius HINTS ${SIRIUS_ROOT}/lib/ ${SIRIUS_ROOT}/lib64/ lib64/)
    # Check if the SIRIUS library was found; if not, report an error.
    if (NOT SIRIUS_LIB)
        message(FATAL_ERROR "SIRIUS build is required but the library was not found.")
    endif()

    get_filename_component(SIRIUS_DIR ${SIRIUS_LIB} DIRECTORY)
    message(STATUS "SIRIUS library found : ${SIRIUS_DIR}")

    # Attempt to find the Fortran module file for SIRIUS.
    find_file(SIRIUS_MOD_FILE NAMES sirius.mod HINTS ${SIRIUS_DIR}/../include/sirius/ ${SIRIUS_DIR}/../include)
    # Check if the module file was found; if not, report an error.
    if (NOT SIRIUS_MOD_FILE)
        message(FATAL_ERROR "SIRIUS Fortran mod file not found")
    endif()

    # Extract the directory of the found Fortran module file.
    get_filename_component(SIRIUS_MOD_DIR ${SIRIUS_MOD_FILE} DIRECTORY)
    message(STATUS "SIRIUS include dirs : ${SIRIUS_MOD_DIR}")
    # Include the directory of the Fortran module for compilation.
    include_directories(${SIRIUS_MOD_DIR})

    # Added SIRIUS definition
    add_compile_definitions(SIRIUS) 

endif()
