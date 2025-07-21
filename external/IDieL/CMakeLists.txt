cmake_minimum_required(VERSION 3.25)

# The version is the atomic number of the current release. Major and minor patches can be set accordingly.
project(IDIEL    VERSION 1.0.0
	         DESCRIPTION "Code for tackling the singularity of the RPA dielectric matrix"
                 LANGUAGES Fortran C CXX)

# Set standard minimum for CXX to 17 (force use of modern compiling framework)
# Dirty trick, alternatively we can make a Fortran 2018 program to check if it is compiled
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#Some options
option(BENCHMARK "Generate the benchmark executables" OFF)
option(INTERFACES "Python bindings" OFF)
option(DOUBLE "Use double precision in real and complex" ON)
option(TESTS "Compile and set the test suit" ON)

# Set precision settings
if(NOT DOUBLE)
    message(STATUS "Compiling with single precision")
    add_compile_definitions(USE_SINGLE_PRECISION)
endif()

# Set default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")

 # Sets the variable in the CACHE. Can view/modify with ccmake
 # FORCE will overwrite if variable is already in the CACHE
 set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Choose the type of build." FORCE)

  # Set the possible values of build type
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

# Set openmp support
include(cmake/openmp.cmake)

# Compiler flags
include(cmake/CompilerFlags.cmake)

# Compiler warnings
include(cmake/CompilerWarnings.cmake)

# Create files
include(cmake/create_files.cmake)

# Intel MKL checks
include(cmake/mkl.cmake)

# FFTW3 library
include(cmake/fftw3.cmake)

# Linear algebra
include(cmake/linalg.cmake)

# Define GNU standard installation directories
include(GNUInstallDirs)

# SPGLIB
include(cmake/spglib.cmake)

# Deviceacc
# Configure the device offloading using deviceacc library
set(DEVICEACC_SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/deviceacc")
include(external/deviceacc/include.cmake)

# For GPU aware compilation we require MPI
if (NOT CPUBACKEND)
    include(cmake/mpi.cmake)
endif()

# Documentation
include(cmake/documentation.cmake)

# Create INSTALL DIRECTORY
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set_property(CACHE CMAKE_INSTALL_PREFIX PROPERTY VALUE "${PROJECT_SOURCE_DIR}/IDieL")
endif()
set(INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
message(STATUS "INSTALL DIRECTORY : ${CMAKE_INSTALL_PREFIX}")


# Set the components source
set(SRC "${SRC_DEVICEACC}")
set(MODDIR "")
add_subdirectory(src)

# Python bindings
if (INTERFACES AND CPUBACKEND)
  file(MAKE_DIRECTORY ${INSTALL_DIR}/interfaces)
  set(PYBIND_CMAKE_DIR "None" CACHE PATH "Path to CMake configuration of pybind")
  message(STATUS "Pybind11 CMake Configuration path : ${PYBIND_CMAKE_DIR}")
  add_subdirectory(interfaces)
elseif(INTERFACES AND NOT CPUBACKEND)
  message(STATUS "Interfaces do not yet support GPU offload: deactivating")
endif()

# Compile test
if(TESTS)
  add_subdirectory(test)
endif()

# If selected compile benchmark programs
# it also downloads all needed data for the benchmark
if(BENCHMARK)
  add_subdirectory(benchmark)
endif(BENCHMARK)

