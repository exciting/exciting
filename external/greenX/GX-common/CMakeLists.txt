# CMakeLists.txt for GreenX Common Modules, required by all other libs in
# the package

add_library(GXCommon STATIC "")

set_target_properties(GXCommon
        PROPERTIES
        VERSION 0.0.1
        SOVERSION 0.0.1
        LIBRARY_OUTPUT_NAME GXCommon
	Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/greenX-common/
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

target_include_directories(GXCommon PUBLIC src/)

target_sources(GXCommon PRIVATE
        src/kinds.f90
        src/constants.f90
        src/unit_conversion.f90
        src/error_handling.f90
        src/lapack_interfaces.f90
        )

