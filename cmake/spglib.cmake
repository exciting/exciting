# This CMake file finds spglib, and in case it 
# is not found it decays to the prepacked version 

option(SPGLIB "Compile with spglib support" ON)
option(USE_INTERNAL_SPGLIB "Use bundled spglib" ON)

if (SPGLIB)
  add_compile_definitions(USE_SPGLIB)

  if (NOT USE_INTERNAL_SPGLIB)
    # Try to find an existing Spglib library
    find_library(SPGLIB_STATIC
      NAMES spg spglib libsymspg symspg
      PATHS ${SPGLIBDIR}
    )

    if (NOT SPGLIB_STATIC)
        message(FATAL_ERROR "External spglib cannot be found. You can use USE_INTERNAL_SPGLIB option to use the bundled version.")
    endif()

    set(spglibInstallDir "")
    add_library(libsymspg STATIC IMPORTED)
    set_target_properties(libsymspg PROPERTIES
      IMPORTED_LOCATION "${SPGLIB_STATIC}"
    )
    message(FATAL_ERROR "NOP")
  else()
    set(spglibInstallDir "${CMAKE_BINARY_DIR}/INTERNAL_spglib_install")
    file(MAKE_DIRECTORY ${spglibInstallDir})
    include(ExternalProject)

    ExternalProject_Add(INTERNAL_SPGLIB
      SOURCE_DIR    "${CMAKE_SOURCE_DIR}/external/spglib-2.6.0/"
      BINARY_DIR    "${CMAKE_BINARY_DIR}/INTERNAL_spglib_build"
      CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER}
                        CFLAGS=-fPIC FCCPP=cpp
                        ${CMAKE_COMMAND} -S <SOURCE_DIR> -B <BINARY_DIR>
                        -DSPGLIB_SHARED_LIBS=OFF
                        -DCMAKE_INSTALL_PREFIX=${spglibInstallDir}
                        -DCMAKE_INSTALL_LIBDIR=lib
                        -DSPGLIB_WITH_TESTS=OFF
      BUILD_COMMAND   $(MAKE)
      INSTALL_COMMAND $(MAKE) install
    )

    add_library(libsymspg STATIC IMPORTED GLOBAL)
    set_target_properties(libsymspg PROPERTIES
      IMPORTED_LOCATION "${spglibInstallDir}/lib/libsymspg.a"
    )
    add_dependencies(libsymspg INTERNAL_SPGLIB)
  endif()
else()
  # Dummy library
  add_library(libsymspg INTERFACE)
endif()
