# exciting Top-level CMakeLists.txt
#
# CMake Variable Descriptions
# PROJECT_BINARY_DIR:  Build directory of the most-recent project() command
# CMAKE_BINARY_DIR: Full path to the top level of the current CMake build tree

# Notice that we need the parser support for Fortran's assumed type, i.e. type(*)
# Older versions (e.g. 3.20) generate a build that at some point asks for a module type within files containing type(*).

# TODO(mrm24): Collapse all exciting (i.e. not in external includes as DEVICEACC options into a single cmake file to be included)
#              Othewise I find very complex to find the options, althoug one can always use CMake to investigate those

cmake_minimum_required(VERSION 3.25)

# The version is the atomic number of the current release. Major and minor patches can be set accordingly.
project(exciting VERSION 11.0.0
                 DESCRIPTION "An all-electron, LAPW, density functional package for ground state calculations and beyond"
                 LANGUAGES Fortran C CXX)

# Set standard minimum for CXX to 17 (force use of modern compiling framework)
# Dirty trick, alternatively we can make a Fortran 2018 program to check if it is compiled
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set folder structure for build directory
set(CMAKE_Fortran_BIN_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_Fortran_LIB_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# Set install directory
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set_property(CACHE CMAKE_INSTALL_PREFIX PROPERTY VALUE "${PROJECT_SOURCE_DIR}/install")
endif()
message(STATUS "Install directory: ${CMAKE_INSTALL_PREFIX}")

# Option set if recompile the external projects
option(RECOMPILE_LIBX_FOX "Recompile libXC and FoX (needed if the compiler changes)" OFF)
if (RECOMPILE_LIBX_FOX)
    set(RECOMPILE_EXT TRUE)
else()
    set(RECOMPILE_EXT FALSE)
endif(RECOMPILE_LIBX_FOX)

# Set the build kind
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build: Debug Release MinSizeRel" FORCE)
endif()

if(NOT CMAKE_BUILD_TYPE MATCHES "Debug|Release")
    message(FATAL_ERROR "Build type ${CMAKE_BUILD_TYPE} is invalid. Valid options are : Debug OR Release. ")
endif()

message(STATUS "Build kind : ${CMAKE_BUILD_TYPE}")

# Set some basic definitions
add_compile_definitions(XS ISO LIBXC)

# Set up MPI
include(cmake/mpi.cmake)

# Set openmp support
include(cmake/openmp.cmake)

# Set the exciting executable name
include(cmake/exciting_name.cmake)

# HDF5 support
include(cmake/hdf5.cmake)

# Intel MKL checks
include(cmake/mkl.cmake)

# FFTW3 library
include(cmake/fftw3.cmake)

# Linear algebra
include(cmake/linalg.cmake)

# Scalapack support
include(cmake/scalapack.cmake)

# SIRIUS support (https://github.com/electronic-structure/SIRIUS/)
include(cmake/sirius.cmake)

# Define GNU standard installation directories
include(GNUInstallDirs)

# CMake module directory
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Defines build type
include(cmake/StandardProjectSettings.cmake)

# Compiler flags
include(cmake/CompilerFlags.cmake)

# Compiler warnings
include(cmake/CompilerWarnings.cmake)

# Configure the device offloading using deviceacc library
set(DEVICEACC_SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/deviceacc")
include(external/deviceacc/include.cmake)

# Create files needed for CMake compilation
include(cmake/create_files.cmake)

# Spglib
include(cmake/spglib.cmake)

# IDieL (it requires all options in previous includes)
include(cmake/idiel.cmake)

# Option controling if we build exciting and needed libraries
# Useful if we only want to build the documentation
option(BUILD_EXCITING "Do we build exciting" ON)

if (BUILD_EXCITING)
    # Compile the needed parts of greenX
    add_subdirectory(external/greenX)
    # External libs, packaged with source
    # LibXC should eventually be switched out for V5
    include(cmake/libxc.cmake)
    include(cmake/fox.cmake)
    add_subdirectory("${PROJECT_SOURCE_DIR}/external/Lbfgsb.3.0")
    add_subdirectory("${PROJECT_SOURCE_DIR}/external/leblaiklib")

    # Build libBZINT
    include(cmake/libbzint.cmake)

    # Location of exciting's CMakeLists
    add_subdirectory(src)

    # Testing
    # Enable ctest
    enable_testing()

    # Testing CMake instructions are found in test folder
    add_subdirectory(test)
endif()

# Documentation
include(cmake/documentation.cmake)


