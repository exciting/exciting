# Find and define build libXC.
#
# Look for version 2.0.0. If not found, build the packaged version with ExternalProject_Add
# Once the libXC API is extended to V5, this can be updated.
#
# Developer Notes
# -----------------
# * ExternalProject_Add argument options should NOT be in quotes.
# * Use add_library t define the library in CMake. find_package necessarily work
#   because there's no guarantee that the library has already been build at
#   CMake configuration time.
# * Useful Stackoverflow ref for linking:
#   https://stackoverflow.com/questions/6351609/cmake-linking-to-library-downloaded-from-externalproject-add

find_package(Libxc 2.0.0 COMPONENTS Fortran)

If (NOT Libxc_FOUND)
  message("-- LibXC not found on the system. Falling back to prepackaged LibXC 2.0.0")

  # Make an installation directory
  set(libXCInstallDir "${CMAKE_SOURCE_DIR}/external/libXC/install_dir")
  file(MAKE_DIRECTORY ${libXCInstallDir})
  message("-- LibXC 2.0.0 will be installed to ${libXCInstallDir}")

  # Build version 2.0.0. packaged with exciting
  include(ExternalProject)
  ExternalProject_Add(INTERAL_LIBXC
    # Location of prepackaged libXC
    SOURCE_DIR    "${CMAKE_SOURCE_DIR}/external/libXC/"
    BUILD_ALWAYS  FALSE
    BUILD_IN_SOURCE 1
    # Based on LibXC's READMe and exciting's make.inc, there should be no need to set
    # C or fortran compiler flags for libXC's production build.
    # But if one wishes to set them, add: CC=<C flags> and/or FCFLAGS=<F90 flags>
    CONFIGURE_COMMAND    ./configure FC=${CMAKE_Fortran_COMPILER} CC=${CMAKE_C_COMPILER} FCCPP=cpp --enable-static=yes --enable-shared=no --prefix=${libXCInstallDir}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )

  add_library(LIBXC STATIC IMPORTED)
  set_target_properties(LIBXC PROPERTIES IMPORTED_LOCATION ${libXCInstallDir}/lib/libxc.a)

  # Headers in include/ are not set as won't be needed for using the static lib with exciting
  # include_directories(${libXCInstallDir}/include)
endif()
