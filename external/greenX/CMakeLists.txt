# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Define GNU standard installation directories
include(GNUInstallDirs)

# Optional multiple precision arithmetic

# Our library directories
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
add_subdirectory(GX-common)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY}/greenX-common/)
add_subdirectory(GX-TimeFrequency)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY}/greenX-MiniMax/)
add_dependencies(LibGXMiniMax GXCommon)
set(GX_MODULES "${CMAKE_Fortran_MODULE_DIRECTORY}/greenX-common/;${CMAKE_Fortran_MODULE_DIRECTORY}/greenX-MiniMax/"  PARENT_SCOPE)
