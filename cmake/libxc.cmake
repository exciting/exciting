# Find and define build libXC.
#
# Look for version 7.0.0 on the system. If it is not found, or if the user asks
# for it, download and build Libxc 7.0.0 from https://gitlab.com/libxc/libxc/
#

# Set the root directory for the libXC installation, allowing for non-standard locations.
# If the user has a custom location for libxc, they can set LIBXC_ROOT.
set(LIBXC_ROOT "None" CACHE STRING "Root directory for non-standard locations for the libxc installation")

# Version, repository and tag used when Libxc is built by exciting
set(LIBXC_DOWNLOAD_VERSION "7.0.0" CACHE STRING "Version of libXC downloaded and built by exciting")
set(LIBXC_DOWNLOAD_REPOSITORY "https://gitlab.com/libxc/libxc.git" CACHE STRING "Git repository the libXC sources are fetched from")

option(USE_INTERNAL_LIBXC "Download and build libXC instead of using a system installation" ON)

if (NOT USE_INTERNAL_LIBXC)
  find_library(LIBXC_LIB NAMES xc HINTS ${LIBXC_ROOT}/lib/)
  find_file(LIBXC_VERSION_FILE NAMES xc_version.h HINTS ${LIBXC_ROOT}/include/)
endif()

if (USE_INTERNAL_LIBXC)
  message("-- LibXC ${LIBXC_DOWNLOAD_VERSION} will be fetched from ${LIBXC_DOWNLOAD_REPOSITORY}")

  # Make an installation directory
  set(libXCInstallDir "${CMAKE_BINARY_DIR}/INTERNAL_libXC_install")
  file(MAKE_DIRECTORY ${libXCInstallDir})
  message("-- LibXC ${LIBXC_DOWNLOAD_VERSION} will be installed to ${libXCInstallDir}")

  if (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    set(FORTRAN_FLAGS_XC "-O0 -ef -g -e Z -dC -s real64 -s integer32 -fPIC -h flex_mp=strict")
  else()
    set(FORTRAN_FLAGS_XC "-cpp ")
  endif()

  # Libxc 7.0.0 requests a CMake policy version below the floor supported by
  # CMake 4, which refuses to configure it unless the minimum is set explicitly
  if (CMAKE_VERSION VERSION_GREATER_EQUAL 4.0)
    set(LIBXC_POLICY_ARGS -DCMAKE_POLICY_VERSION_MINIMUM=3.5)
  endif()

  # Download and build libXC
  include(ExternalProject)
  ExternalProject_Add(INTERNAL_LIBXC
    # Upstream location of libXC
    GIT_REPOSITORY ${LIBXC_DOWNLOAD_REPOSITORY}
    GIT_TAG        ${LIBXC_DOWNLOAD_VERSION}
    GIT_SHALLOW    TRUE
    GIT_PROGRESS   TRUE
    # Do not contact the server again once the sources have been fetched
    UPDATE_COMMAND ""
    SOURCE_DIR    "${CMAKE_BINARY_DIR}/INTERNAL_libXC_source"
    BUILD_ALWAYS   ${RECOMPILE_EXT}
    BINARY_DIR    "${CMAKE_BINARY_DIR}/INTERNAL_libXC_build"
    CONFIGURE_COMMAND FC=${CMAKE_Fortran_COMPILER} CC=${CMAKE_C_COMPILER}
                      CFLAGS=-fPIC FCCPP=cpp 
                      ${CMAKE_COMMAND} -H<SOURCE_DIR> -B<BINARY_DIR>
                      -DCMAKE_INSTALL_PREFIX=${libXCInstallDir}
                      -DBUILD_TESTING=OFF
                      ${LIBXC_POLICY_ARGS}
		      -DCMAKE_INSTALL_LIBDIR=lib
                      -DENABLE_FORTRAN=ON
		      -DCMAKE_Fortran_FLAGS=${FORTRAN_FLAGS_XC}
    BUILD_COMMAND   $(MAKE) clean && $(MAKE)
    INSTALL_COMMAND $(MAKE) install
  )

  include_directories(${libXCInstallDir}/include)
  set(LIBXC_LIBDIR "${libXCInstallDir}/lib")

  add_library(XC INTERFACE)
  target_link_libraries(XC INTERFACE ${LIBXC_LIBDIR}/libxcf03.a ${LIBXC_LIBDIR}/libxc.a)
  add_dependencies(XC INTERNAL_LIBXC)
  add_compile_definitions(LIBXC_HAS_FUNC_MOD)

elseif(NOT LIBXC_LIB AND NOT USE_INTERNAL_LIBXC)
  message(FATAL_ERROR "External Libxc was required but not found")

else()
  # Read the contents of the file
  file(READ "${LIBXC_VERSION_FILE}" LIBXC_HEADER_CONTENTS)
  # Extract XC_VERSION
  string(REGEX MATCH "#define XC_VERSION \"([0-9]+\\.[0-9]+\\.[0-9]+)\"" _match "${LIBXC_HEADER_CONTENTS}")
  set(LIBXC_VERSION "${CMAKE_MATCH_1}")
  # Extract major, minor, and micro versions
  string(REGEX MATCH "#define XC_MAJOR_VERSION ([0-9]+)" _match_major "${LIBXC_HEADER_CONTENTS}")
  set(LIBXC_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "#define XC_MINOR_VERSION ([0-9]+)" _match_minor "${LIBXC_HEADER_CONTENTS}")
  set(LIBXC_MINOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "#define XC_MICRO_VERSION ([0-9]+)" _match_micro "${LIBXC_HEADER_CONTENTS}")
  set(LIBXC_MICRO_VERSION "${CMAKE_MATCH_1}")

  # Display the version info (optional)
  message(STATUS "LibXC Version: ${LIBXC_VERSION}")
  if (LIBXC_MAJOR_VERSION LESS 5)
    message(FATAL_ERROR "exciting does only support libxc version 5.0.0 or higher.")
  endif()
  get_filename_component(Libxc_libdir ${LIBXC_LIB} DIRECTORY)
  find_library(xcf03 NAMES xcf03 HINTS ${LIBXC_ROOT}/lib/)
  message(STATUS "Libxc external libs: ${LIBXC_LIB}:${xcf03}")
  get_filename_component(Libxc_includes ${LIBXC_VERSION_FILE} DIRECTORY)
  include_directories(${Libxc_includes})

  # We check that the installation contain Fortran wrappers
  file(GLOB FILES_IN_DIR "${Libxc_includes}/*")
  foreach(FILE ${FILES_IN_DIR})
    string(TOLOWER "${FILE}" LOWER_FILE)
    if(LOWER_FILE MATCHES "xc_f03_lib_m\\.mod$")
        set(FILE_FOUND TRUE)
        break()
    endif()
  endforeach()
  if(NOT FILE_FOUND)
    message(FATAL_ERROR "Libxc does not contain Fortran wrappers.")
  endif()

  add_library(XC INTERFACE)
  target_link_libraries(XC INTERFACE ${xcf03} ${LIBXC_LIB})
  add_dependencies(XC INTERNAL_LIBXC)
  
  # For libXC 7.0.0 there was an API change. Internally we control those changes using pragmas.
  if (LIBXC_MAJOR_VERSION GREATER_EQUAL 7)
    add_compile_definitions(LIBXC_HAS_FUNC_MOD)
  endif()
endif()

