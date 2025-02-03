cmake_minimum_required(VERSION 3.15)

project(LEBLAIK_PROJ VERSION 1.0.0 LANGUAGES Fortran)

# Create the library target
add_library(leblaik STATIC "")

# General compiler options for the target
if (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    target_compile_options(leblaik PRIVATE -O2 -cpp -fPIC -s real64 -s integer32)
else()
    target_compile_options(leblaik PRIVATE -O3 -cpp -fPIC)
endif()

# Set target properties
set_target_properties(leblaik
    PROPERTIES
    VERSION 1.0.0
    LINKER_LANGUAGE Fortran
    LIBRARY_OUTPUT_NAME leblaik
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
)

# Specify source files
set(free_form_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/leblaik.f90
)

set(fixed_form_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/Lebedev-Laikov.F
)

# Add source files to the target
target_sources(leblaik PRIVATE
    ${free_form_sources}
    ${fixed_form_sources}
)

# Set compile options for fixed-form Fortran source files
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        foreach(source_file ${fixed_form_sources})
                set_source_files_properties(${source_file} PROPERTIES
                        COMPILE_OPTIONS "-fixed"
                )
        endforeach()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang")
        foreach(source_file ${fixed_form_sources})
                set_source_files_properties(${source_file} PROPERTIES
                        COMPILE_OPTIONS "-ffixed-form"
                )
        endforeach()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
        foreach(source_file ${fixed_form_sources})
                set_source_files_properties(${source_file} PROPERTIES
                        COMPILE_OPTIONS "-ffixed; -s integer32; -s real64; -h flex_mp=strict"
                )
        endforeach()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        foreach(source_file ${fixed_form_sources})
                set_source_files_properties(${source_file} PROPERTIES
                        COMPILE_OPTIONS "-std=legacy"
                )
        endforeach()
else()
        message(FATAL_ERROR "Unrecognized compiler: only Intel, GNU, Cray, Flang (new) compilers are supported")
endif()

