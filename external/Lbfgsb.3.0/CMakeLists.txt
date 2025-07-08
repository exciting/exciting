cmake_minimum_required(VERSION 3.15)

project(LBFGSB_PROJ VERSION 3.0
        LANGUAGES Fortran)

add_library(lbfgsb STATIC "")

if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        target_compile_options(lbfgsb PRIVATE -O3 -fPIC -fixed)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang")
        target_compile_options(lbfgsb PRIVATE -O3 -fPIC -ffixed-form)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
        target_compile_options(lbfgsb PRIVATE -O3 -fPIC -ffixed -s integer32 -s real64)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        target_compile_options(lbfgsb PRIVATE -O3 -std=legacy -fPIC)
else()
        message(FATAL_ERROR "Unrecognized compiler: only Intel, GNU, Cray, Flang (new) compilers are supported")
endif()

set_target_properties(lbfgsb
        PROPERTIES
        VERSION 3.0
    	LINKER_LANGUAGE Fortran
        LIBRARY_OUTPUT_NAME lbfgsb
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

target_sources(lbfgsb PRIVATE
	timer.f
	linpack.f
	lbfgsb.f
	)
