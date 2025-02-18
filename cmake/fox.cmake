# Build FoX library

set(FoXInstallDir "${CMAKE_BINARY_DIR}/INTERNAL_FoX_install")
file(MAKE_DIRECTORY ${FoXInstallDir})

message(STATUS "Using prepackaged FoX library for XML parsing.")
message(STATUS "FoX (version 2012) will be installed to ${FoXInstallDir}")

# Set compiler-specific flags
if (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    set(FORTRAN_FLAGS_FOX "-O0 -ef -hipa0 -g -e Z -dC -s real64 -s integer32 -fPIC -h flex_mp=strict")
else()
    set(FORTRAN_FLAGS_FOX "-fPIC")
endif()

# Configure and build the FoX library
include(ExternalProject)
ExternalProject_Add(INTERNAL_FOX
  # Location of prepackaged libXC
  SOURCE_DIR    "${CMAKE_SOURCE_DIR}/external/FoX/"
  BUILD_ALWAYS   ${RECOMPILE_EXT}
  BINARY_DIR    "${CMAKE_BINARY_DIR}/INTERNAL_FoX_build"
  CONFIGURE_COMMAND FC=${CMAKE_Fortran_COMPILER}
                    ${CMAKE_COMMAND} -H<SOURCE_DIR> -B<BINARY_DIR>
                    -DCMAKE_INSTALL_PREFIX=${FoXInstallDir}
                    -DFoX_ENABLE_EXAMPLES=OFF
                    -DCMAKE_INSTALL_LIBDIR=lib
                    -DFoX_ENABLE_WKML=OFF
                    -DCMAKE_Fortran_FLAGS=${FORTRAN_FLAGS_XC}
  BUILD_COMMAND   $(MAKE) clean && $(MAKE)
  INSTALL_COMMAND $(MAKE) install
)

# Add FoX include and library directories
include_directories(${FoXInstallDir}/include)
link_directories(${FoXInstallDir}/lib)

# Helper macro to create imported FoX libraries
macro(add_FoX_library target_name lib_name)
    add_library(${target_name} STATIC IMPORTED)
    set_property(TARGET ${target_name} PROPERTY IMPORTED_LOCATION ${FoXInstallDir}/lib/lib${lib_name}.a)
    add_dependencies(${target_name} INTERNAL_FOX)
endmacro()

# Define imported FoX libraries
add_FoX_library(FoX_dom     FoX_dom)
add_FoX_library(FoX_sax     FoX_sax)
add_FoX_library(FoX_utils   FoX_utils)
add_FoX_library(FoX_wxml    FoX_wxml)
add_FoX_library(FoX_wcml    FoX_wcml)
add_FoX_library(FoX_fsys    FoX_fsys)
add_FoX_library(FoX_common  FoX_common)
