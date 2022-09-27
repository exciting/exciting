cmake_minimum_required(VERSION 3.15)

project(LBFGSB VERSION 3.0.0
        LANGUAGES Fortran)

add_library(liblbfgsb STATIC "")

target_compile_options(liblbfgsb PRIVATE -03 -std=legacy)

set_target_properties(liblbfgsb
        PROPERTIES
	    LINKER_LANGUAGE Fortran
        VERSION 3.0
        LIBRARY_OUTPUT_NAME lbfgsb
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

target_sources(liblbfgsb PRIVATE
	timer.f
	linpack.f
	lbfgsb.f
	)
