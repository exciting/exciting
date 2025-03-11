# Build an exciting library with it
# Compile the exciting executable and link to the exciting library 

###############################################
# Get all source objects for exciting library #
###############################################

# Define the current folder
set(EXCITINGLIB_SRC_FOLDERS ${CMAKE_CURRENT_SOURCE_DIR})
set(EXCITINGLIB_SRC_FILES "")
set(EXCITINGLIB_SRC_FILES_FIXED "")

# Specify folders to ignore
set(FOLDER_SPACEGROUP "${CMAKE_CURRENT_SOURCE_DIR}/spacegroup")
set(FOLDER_LIBBZINT   "${CMAKE_CURRENT_SOURCE_DIR}/src_libbzint")

# Specify files not inside ignored folders to ignore
# TODO: Remove or check those files
set(EXCLUDE_SRC_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/species/species.f90"
    # Add more files to exclude as needed
)

if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel" AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    message(STATUS "Classic Intel: using old GW files due to c_f_pointer bug.")
    list(APPEND EXCLUDE_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src_gw/calcminm2.f90")
    list(APPEND EXCLUDE_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src_gw/calcepsilon.f90")
    list(APPEND EXCLUDE_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src_gw/calcmwm.f90")
    list(APPEND EXCLUDE_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src_gw/expand_products.f90")
else()
    list(APPEND EXCLUDE_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src_gw/calcminm2_ifort.f90")
    list(APPEND EXCLUDE_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src_gw/calcepsilon_ifort.f90")
    list(APPEND EXCLUDE_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src_gw/calcmwm_ifort.f90")
    list(APPEND EXCLUDE_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src_gw/expand_products_ifort.f90")
endif()

# We need special treatment of files containing INQUIRE so that it compiles using the Cray compiler 
# Developers, please add here anycase of INQUIRE using not the default integer
# TODO: remove when Cray compiler allows for arbitrary precision recl 
set(EXCITING_INQUIRE_FILES "")
if (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")

    list(APPEND EXCITING_INQUIRE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src_xs/putgeteps0.f90")
    list(APPEND EXCITING_INQUIRE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src_xs/src_rttddft/rttddft_io_serial.f90")
    list(APPEND EXCITING_INQUIRE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src_gw/calcpmatgw.f90")
    list(APPEND EXCITING_INQUIRE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src_hybrids/putvxnl.f90")

    # We create an INQUIRE_FLAGS
    set(INQUIRE_FLAGS "-s integer64")
    set_source_files_properties(${EXCITING_INQUIRE_FILES} PROPERTIES COMPILE_FLAGS ${INQUIRE_FLAGS})
endif()

# Iterate over files and purge main
foreach(folder ${EXCITINGLIB_SRC_FOLDERS})
    # Search for F90 files in the current folder
    file(GLOB_RECURSE F90_FILES_IN_FOLDER "${folder}/*.f90" "${folder}/*.F90")

    # Exclude folder as this
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "^${FOLDER_SPACEGROUP}/.*$")
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "^${FOLDER_LIBBZINT}/.*$")

    # Iterate over files and exclude those in EXCLUDE_FILES list
    foreach(file ${EXCLUDE_SRC_FILES})
        list(FIND F90_FILES_IN_FOLDER "${file}" EXCLUDE_FILE_INDEX)
        if(EXCLUDE_FILE_INDEX GREATER -1)
            list(REMOVE_AT F90_FILES_IN_FOLDER ${EXCLUDE_FILE_INDEX})
        endif()
    endforeach()

    # Iterate over files and exclude those in INQUIRE_FILES list
    foreach(file ${EXCITING_INQUIRE_FILES})
        list(FIND F90_FILES_IN_FOLDER "${file}" EXCLUDE_FILE_INDEX)
        if(EXCLUDE_FILE_INDEX GREATER -1)
            list(REMOVE_AT F90_FILES_IN_FOLDER ${EXCLUDE_FILE_INDEX})
        endif()
    endforeach()

    # Exclude main.f90 and main.F90
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "(^|/)main\\.f90$")
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "(^|/)main\\.F90$")

    # Append the filtered files to the list
    list(APPEND EXCITINGLIB_SRC_FILES ${F90_FILES_IN_FOLDER})
endforeach()

# Now process the old fixed format Fortran files
foreach(folder ${EXCITINGLIB_SRC_FOLDERS})
    # Search for F90 files in the current folder
    file(GLOB_RECURSE F90_FILES_IN_FOLDER "${folder}/*.f")

    # Exclude folder as this
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "^${FOLDER_SPACEGROUP}/.*$")
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "^${FOLDER_LIBBZINT}/.*$")

    # Iterate over files and exclude those in EXCLUDE_FILES list
    foreach(file ${EXCLUDE_SRC_FILES})
        list(FIND F90_FILES_IN_FOLDER "${file}" EXCLUDE_FILE_INDEX)
        if(EXCLUDE_FILE_INDEX GREATER -1)
            list(REMOVE_AT F90_FILES_IN_FOLDER ${EXCLUDE_FILE_INDEX})
        endif()
    endforeach()

    # Iterate over files and exclude those in INQUIRE_FILES list
    foreach(file ${EXCITING_INQUIRE_FILES})
        list(FIND F90_FILES_IN_FOLDER "${file}" EXCLUDE_FILE_INDEX)
        if(EXCLUDE_FILE_INDEX GREATER -1)
            list(REMOVE_AT F90_FILES_IN_FOLDER ${EXCLUDE_FILE_INDEX})
        endif()
    endforeach()


    # Exclude main.f90 and main.F90
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "(^|/)main\\.f90$")
    list(FILTER F90_FILES_IN_FOLDER EXCLUDE REGEX "(^|/)main\\.F90$")

    # Append the filtered files to the list
    list(APPEND EXCITINGLIB_SRC_FILES_FIXED ${F90_FILES_IN_FOLDER})
endforeach()

# IntelLLVM (version 2025.0.0) compiler fails for some reason to automatically set 
# fixed form for *.f files.
if (CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
        foreach(source_file ${EXCITINGLIB_SRC_FILES_FIXED})
                set_source_files_properties(${source_file} PROPERTIES
                        COMPILE_OPTIONS "-fixed"
                )
        endforeach()
endif()

# For Cray compiler get further info about the compilation of complex
# GPU files. It only affects GPU-aware compilations.
if (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    set(MARK_FILES "")
    list(APPEND MARK_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src_gw/calcminm2.f90")

    # We create an MARK_FLAGS
    set(MARK_FLAGS "-hlist=a ")
    set_source_files_properties(${MARK_FILES} PROPERTIES COMPILE_FLAGS ${MARK_FLAGS})

endif()

# Developers: uncomment that to check if the source files that you
# want in the project to be included
# message("Source files:")
# foreach(file ${EXCITINGLIB_SRC_FILES})
#    message("${file}")
# endforeach()
# message(FATAL_ERROR "EXIT")

#####################################
# Compile and link exciting library #
#####################################

add_library(excitingLib_objects OBJECT ${SRC_DEVICEACC} ${EXCITINGLIB_SRC_FILES} ${EXCITINGLIB_SRC_FILES_FIXED} ${EXCITING_INQUIRE_FILES})
                                    
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
# Check exciting has all dependencies compiled
add_dependencies(excitingLib_objects XC FoX_dom FoX_sax FoX_utils FoX_wcml FoX_wxml FoX_fsys FoX_common)

# Here create exciting library
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    add_library(excitingLib SHARED $<TARGET_OBJECTS:excitingLib_objects>)
else()
    add_library(excitingLib STATIC $<TARGET_OBJECTS:excitingLib_objects>)
endif()
target_link_directories(excitingLib PUBLIC ${CMAKE_Fortran_LIB_DIRECTORY})

#Create and interface for FoX
add_library(FoX INTERFACE)
target_link_libraries(FoX INTERFACE FoX_dom FoX_sax FoX_wxml FoX_common FoX_utils FoX_fsys)

#Link exciting library with other libs
target_link_libraries(excitingLib PUBLIC leblaik lbfgsb XC GXCommon LibGXMiniMax bzint
	                          FoX ${SCALAPACK_LIB} ${SIRIUS_LIB}
                                  ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${HDF5_LIBS}
                                  ${FFTW3_LIBRARIES} ${MPI_LIBS} ${LINKS_DEVICEACC})

# Set the linker to Fortran, otherwise C++ might be used in some cases
set_property(TARGET excitingLib PROPERTY LINKER_LANGUAGE Fortran)

#############################
# Build exciting executable #
#############################

add_executable(${exciting_name} "${CMAKE_CURRENT_SOURCE_DIR}/mainxml/main.f90")

# Linking 
target_link_directories(${exciting_name} PUBLIC ${CMAKE_Fortran_LIB_DIRECTORY})
target_link_libraries(${exciting_name} PUBLIC excitingLib)

# Sometimes GNU do not resolve the dependency on C++
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_link_libraries(${exciting_name} PRIVATE stdc++)
endif()


# For Intel and GNU indicate that exciting has the objects available for LTO
# In the rest of compilers this is done with compiler flags or for Cray is fully deactivated
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" OR CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    set_property(TARGET ${exciting_name} PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
endif()
# Set the linker to Fortran, otherwise C++ linker is used
set_property(TARGET ${exciting_name} PROPERTY LINKER_LANGUAGE Fortran)

# Ensure libraries are found

if(MAGMA_DIR)
    set_target_properties(${exciting_name} PROPERTIES
                INSTALL_RPATH "${MAGMA_DIR}/lib/libmagma.so"
                INSTALL_RPATH_USE_LINK_PATH TRUE)
else()
    set_target_properties(${exciting_name} PROPERTIES
                INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

# Install modules in case one wants to use excitinglib to build anything else
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/
        DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_EXECUTE
                            WORLD_READ WORLD_EXECUTE
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE
                         GROUP_READ
                         WORLD_READ
)

# Install exciting library
set_target_properties(excitingLib PROPERTIES
    INSTALL_RPATH "${MAGMA_DIR}/lib/"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

install(TARGETS excitingLib
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

# Set rpath for the exciting
set_target_properties(${exciting_name} PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:${MAGMA_DIR}/lib/"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

# Install exciting
install(TARGETS ${exciting_name}
        DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)

# Create the usual exciting wrapper for simplified calls
install(CODE "execute_process(COMMAND ln -sf ${CMAKE_INSTALL_PREFIX}/bin/${exciting_name} ${CMAKE_INSTALL_PREFIX}/bin/exciting)")
