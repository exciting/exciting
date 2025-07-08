# Local components
set(SRC_MEMORY ${CMAKE_CURRENT_SOURCE_DIR}/common/memory_device.f90)

if(NOT CPUBACKEND)
    set(SRC_MEMORY ${SRC_MEMORY}
                  ${CMAKE_CURRENT_SOURCE_DIR}/device/s_memory_device.f90)
else()
    set(SRC_MEMORY ${SRC_MEMORY}
        ${CMAKE_CURRENT_SOURCE_DIR}/host/s_memory_device.f90)
endif()

# Push back to SRC
set(SRC ${SRC} ${SRC_MEMORY} PARENT_SCOPE)
