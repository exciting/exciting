cmake_minimum_required(VERSION 3.13)

project(vdWCorrections VERSION 1.0.0
                       DESCRIPTION "An efficient implementation of van der Waals density 
                                    functionals based on a Monte-Carlo integration technique"
                       LANGUAGES Fortran)

set(CMAKE_Fortran_BIN_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_Fortran_LIB_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# Set library location
set(DEFAULT_CUBA_PATH "${CMAKE_SOURCE_DIR}/Cuba-3.2")
set(CUBA_PATH "${DEFAULT_CUBA_PATH}" CACHE STRING "Location of Cuba-3.2 library")

# Find external library
find_library(CUBA NAME "libcuba" "cuba" HINTS "${CUBA_PATH}")

# If external library is not built, build it
if(NOT CUBA)
  include(ExternalProject)
  ExternalProject_Add(CUBA_EXTERNAL
                    SOURCE_DIR ${DEFAULT_CUBA_PATH}
                    CONFIGURE_COMMAND ./configure --with-maxdim=6 --with-maxcomp=1 --prefix=${DEFAULT_CUBA_PATH}
		    BUILD_COMMAND make
		    BUILD_IN_SOURCE TRUE
                    )	
  set(CUBA "${DEFAULT_CUBA_PATH}/libcuba.a")
endif()

# Main executable 
add_executable(noloco)

set_target_properties(noloco
        PROPERTIES
        RUNTIME_OUTPUT_NAME noloco.exe
        VERSION 1.0.0
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_Fortran_BIN_DIRECTORY}
        )

target_sources(noloco
  PRIVATE
  src/blend103.f90
  src/integrand.f90
  src/main.f90
  src/modules.f90
  src/r3minv.f90
  src/read_densities.f90
  src/read_gradients.f90
  src/read_input.f90
  src/read_phi.f90
  src/vdWDF_LL.f90
  src/vdWDF_VV.f90
        )

target_link_libraries(noloco ${CUBA})
