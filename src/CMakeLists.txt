#####################################
# Compile and link exciting library #
#####################################

add_library(excitingLib_objects OBJECT ${SRC_DEVICEACC} ${EXCITING_SOURCE_FILES})
                                    
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

target_link_libraries(excitingLib PUBLIC leblaik lbfgsb XC GXCommon LibGXMiniMax bzint
	                          FoX ${SCALAPACK_LIB} ${SIRIUS_LIB}
                                  ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${HDF5_LIBS}
                                  ${FFTW3_LIBRARIES} ${MPI_LIBS})

target_link_libraries(excitingLib PUBLIC ${LINKS_DEVICEACC})
# Set the linker to Fortran, otherwise C++ linker is used
set_property(TARGET excitingLib PROPERTY LINKER_LANGUAGE Fortran)

#############################
# Build exciting executable #
#############################

add_executable(${exciting_name} ${EXCITING_MAIN})
target_link_directories(${exciting_name} PUBLIC ${CMAKE_Fortran_LIB_DIRECTORY})
target_link_libraries(${exciting_name} PUBLIC excitingLib)

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



