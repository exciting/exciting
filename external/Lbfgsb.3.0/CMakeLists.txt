cmake_minimum_required(VERSION 3.15)

project(LBFGSB VERSION 3.0
        LANGUAGES Fortran)

add_library(liblbfgsb STATIC "")

target_compile_options(liblbfgsb PRIVATE -O3 -std=legacy)

set_target_properties(liblbfgsb
        PROPERTIES
        VERSION 3.0
    	LINKER_LANGUAGE Fortran
        LIBRARY_OUTPUT_NAME lbfgsb
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

target_sources(liblbfgsb PRIVATE
	timer.f
	linpack.f
	lbfgsb.f
	)
