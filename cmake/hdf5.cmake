# hdf5.cmake
#
# This CMake script enables or disables HDF5 support based on the user-defined
# option HDF5. If HDF5 support is enabled, it finds the required HDF5 
# packages and sets the necessary compilation and linking parameters.
# Additionally checks for parallel IO support. 

# Option to enable or disable HDF5 support
option(HDF5 "Enables HDF5 support" OFF)

# Check if HDF5 support is enabled
if(HDF5)
    # If enabled, display a status message
    message(STATUS "HDF5 support : ON")

    # Set looking for parallel HDF5 version
    set(HDF5_PREFER_PARALLEL ON)

    # Find the HDF5 package, requiring components for C, C++, and Fortran
    find_package(HDF5 COMPONENTS C Fortran REQUIRED)

    # If found but does not support parallel IO end the compilation
    if( NOT HDF5_IS_PARALLEL )
        message(WARNING "For a performant binary exciting requires HDF5 with parallel IO support")
    endif()

    # Define a compilation flag for HDF5
    add_compile_definitions(_HDF5_)

    # Add the HDF5 definitions to the compiler options
    add_definitions(${HDF5_DEFINITIONS})

    # Include directories for HDF5 header files
    include_directories(${HDF5_INCLUDE_DIRS})

    # Link directories for HDF5 libraries
    link_directories(${HDF5_LIBRARY_DIRS})

    # Set the libraries to link against, combining C, C++, and Fortran libraries
    set(HDF5_LIBS "${HDF5_C_LIBRARIES};${HDF5_CXX_LIBRARIES};${HDF5_Fortran_LIBRARIES}")
else()
    # If HDF5 support is not enabled, display a status message
    message(STATUS "HDF5 support : OFF")
endif()
