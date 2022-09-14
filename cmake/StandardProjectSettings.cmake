
# Check compiler is supported
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  if (NOT ${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER "6.9.9" )
    message(FATAL_ERROR "The minimum supported version of GCC is 7.0.0")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  #TODO(Alex) Add minimum version for Intel
endif()

# Set default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(
    STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")

 # Sets the variable in the CACHE. Can view/modify with ccmake
 # FORCE will overwrite if variable is already in the CACHE
 set(CMAKE_BUILD_TYPE
      RelWithDebInfo
      CACHE STRING "Choose the type of build." FORCE)

  # Set the possible values of build type
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel"
  "RelWithDebInfo")
endif()
