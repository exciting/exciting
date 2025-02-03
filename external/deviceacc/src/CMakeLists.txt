# Add the source in each subfolder
add_subdirectory(control)
add_subdirectory(linalg)
add_subdirectory(fft)
add_subdirectory(memory)
add_subdirectory(macros)

# Compile
message(STATUS "FILES: ${SRC}")
set(MOD_DIR "${CMAKE_BINARY_DIR}/include")
file(MAKE_DIRECTORY ${MOD_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${MOD_DIR})
if(CMAKE_Fortran_COMPILER MATCHES "ifx")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -module ${MOD_DIR}")
endif()
message(STATUS "MODULES dir : ${MOD_DIR}")
add_library(deviceacc_objects OBJECT ${SRC})

if(NOT CPUBACKEND)
    # Add symbols to the static library so we only need to link against it and not its dependencies
    if (NVIDIA)
	add_library(deviceacc STATIC $<TARGET_OBJECTS:deviceacc_objects>)
        set(LINKS ${devacc_link_libs} ${MAGMA_DIR}/lib/libmagma.so CUDA::cufft CUDA::cudart)
    elseif(AMD)
        add_library(deviceacc STATIC $<TARGET_OBJECTS:deviceacc_objects>)
        set(LINKS ${devacc_link_libs} ${MAGMA_DIR}/lib/libmagma.so ${rocfftlib})
    elseif(INTEL)
        add_library(deviceacc STATIC $<TARGET_OBJECTS:deviceacc_objects>)
        set(LINKS ${devacc_link_libs})
    else()
        message(FATAL_ERROR "Should not be here")
    endif()

    target_link_libraries(deviceacc PUBLIC ${LINKS})

    if (NVIDIA OR AMD)
        set_target_properties(deviceacc PROPERTIES
            INSTALL_RPATH "${MAGMA_DIR}/lib/libmagma.so"
            INSTALL_RPATH_USE_LINK_PATH TRUE)
    endif()

else()
    add_library(deviceacc STATIC $<TARGET_OBJECTS:deviceacc_objects>)
    set(LINKS ${devacc_link_libs})
    target_link_libraries(deviceacc PUBLIC ${LINKS})
endif()

set_property(TARGET deviceacc_objects PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
set_property(TARGET deviceacc PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

install(DIRECTORY ${MOD_DIR}
        DESTINATION ${INSTALL_DIR}
        DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_EXECUTE
                            WORLD_READ WORLD_EXECUTE
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE
                         GROUP_READ
                         WORLD_READ
        )
install(FILES ${MACROS} DESTINATION ${INSTALL_DIR}/include/)
install(TARGETS deviceacc DESTINATION ${INSTALL_DIR}/lib)
