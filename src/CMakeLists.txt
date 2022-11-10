# CMakeLists.txt for exciting

add_library(excitingLib "")

set_target_properties(excitingLib
        PROPERTIES
        VERSION 8.0.0
        SOVERSION 8.0.0
        LIBRARY_OUTPUT_NAME exciting
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

# TODO Add an if statement here that appends ${COMPILER_WARNINGS} flags
target_compile_options(excitingLib PRIVATE ${STD_FFLAGS})

# TODO Add all source files here, or add CMakeLists.txt in their respective subfolders
target_sources(excitingLib PRIVATE
    constants/precision.f90
    constants/constants.f90)

target_link_libraries(excitingLib liblbfgsb libleblaik ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# TODO Link excitingLib to exciting executable

# TODO Set installation
# Install library
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in the top-level
# CMakeLists.txt
#install(TARGETS excitingLib
#        ARCHIVE DESTINATION lib
#        LIBRARY DESTINATION lib)

# Install modules
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in the top-level
## CMakeLists.txt
#install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
#        DESTINATION include)
