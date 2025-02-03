# **************************************************************************************************
#  Copyright (C) 2020-2023 GreenX library
#  This file is distributed under the terms of the APACHE2 License.
#
# **************************************************************************************************

# CMakeLists for TimeFrequency Library (MiniMax)

# GreenXMiniMax Library
add_library(LibGXMiniMax STATIC "")

# Set properties of the library
set_target_properties(LibGXMiniMax
        PROPERTIES
        VERSION 0.0.1
        SOVERSION 0.0.1
        LIBRARY_OUTPUT_NAME gx_minimax
	Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/greenX-MiniMax/
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

# Add include directories to the LibGXMiniMax
target_include_directories(LibGXMiniMax PUBLIC src/ api/ utilities/)

# Define source that comprise LibGXMiniMax
target_sources(LibGXMiniMax PRIVATE
        src/gx_common.h
	src/minimax_utils.F90
        src/minimax_tau.F90
        src/minimax_omega.F90
        src/minimax_grids.F90
        api/api_utilities.f90
        api/gx_minimax.f90
        )

# Link external libraries and lower-level GX lib dependency
target_link_libraries(LibGXMiniMax GXCommon ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})


