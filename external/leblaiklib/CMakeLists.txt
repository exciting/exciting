cmake_minimum_required(VERSION 3.15)

project(LEBLAIK VERSION 1.0.0
        LANGUAGES Fortran)

add_library(libleblaik STATIC "")

target_compile_options(libleblaik PRIVATE -O3 -std=legacy)

set_target_properties(libleblaik
        PROPERTIES
        VERSION 1.0.0
        LINKER_LANGUAGE Fortran
        LIBRARY_OUTPUT_NAME leblaik
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

target_sources(libleblaik PRIVATE
	leblaik.f90
	Lebedev-Laikov.F
        )
