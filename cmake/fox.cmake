# Build FoX library.

# Make an installation directory
# Note, by default the modules and libs are copied to FOXROOT/objs finclude/ and lib/
# FOXROOT/objs is not the installation path though
set(FoXInstallDir "${CMAKE_SOURCE_DIR}/external/FoX/install")
file(MAKE_DIRECTORY ${FoXInstallDir})

message("-- Using prepackaged FoX library for XML parsing.")
message("-- FoX 2009 will be installed to ${FoXInstallDir}")

# Configure and build
include(ExternalProject)
ExternalProject_Add(INTERAL_FOX
   # Location of prepackaged FoX
   SOURCE_DIR    "${CMAKE_SOURCE_DIR}/external/FoX/"
   BUILD_ALWAYS  FALSE
   BUILD_IN_SOURCE 1
   CONFIGURE_COMMAND  ./configure FC=${CMAKE_Fortran_COMPILER} --prefix=${FoXInstallDir}
   BUILD_COMMAND make
   INSTALL_COMMAND make install
)

# Approach to generating a FOX lib. Note, there are several in ${FoXInstallDir}/lib/
# This may require several targets.

# add_library(FOX STATIC IMPORTED)
# set_target_properties(FOX PROPERTIES IMPORTED_LOCATION ${FoXInstallDir}/lib/)
# include_directories(${libXCInstallDir}/finclude)
