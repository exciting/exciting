#######################################
# Create objects or links using CMake #
#######################################

###################################################################################
# In Cray 18 there is an incompatibility with what is generated and CMake expects #
# for Fortran submodules                                                          #
###################################################################################

if (CMAKE_VERSION VERSION_LESS_EQUAL "4.0.0" AND CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    string(REGEX MATCH "^([0-9]+)" COMPILER_MAJOR_VERSION "${CMAKE_Fortran_COMPILER_VERSION}")
    message(STATUS "Correcting ${CMAKE_Fortran_COMPILER_VERSION}")
    if (COMPILER_MAJOR_VERSION AND COMPILER_MAJOR_VERSION VERSION_GREATER_EQUAL "18")
	# Create the modules folder
	file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/include/)
	# Create the symlink
	execute_process(
		COMMAND ln -s ${CMAKE_BINARY_DIR}/include/m_memory_device.SM_MEMORY_DEVICE.smod
	                      ${CMAKE_BINARY_DIR}/include/sm_memory_device.mod
	
        )
        execute_process(
                COMMAND ln -s ${CMAKE_BINARY_DIR}/include/idiel.idIEL_2D.smod
                              ${CMAKE_BINARY_DIR}/include/idiel_2d.mod

        )
        execute_process(
                COMMAND ln -s ${CMAKE_BINARY_DIR}/include/idiel.idIEL_3D.smod
                              ${CMAKE_BINARY_DIR}/include/idiel_3d.mod

        )
        execute_process(
		COMMAND ln -s ${CMAKE_BINARY_DIR}/include/idiel.idiel_COMMON.smod
                              ${CMAKE_BINARY_DIR}/include/idiel_common.mod

        )
    endif()
endif()

