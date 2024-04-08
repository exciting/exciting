# Local components
if(NOT CPUBACKEND)
    set(SRC_CONTROL ${CMAKE_CURRENT_SOURCE_DIR}/device/device_host_register.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/common/device_host_register.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/device/device_world_t.f90)
else()
    set(SRC_CONTROL ${CMAKE_CURRENT_SOURCE_DIR}/host/device_host_register.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/common/device_host_register.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/host/device_world_t.f90)
endif()
# Push back to SRC
set(SRC ${SRC} ${SRC_CONTROL} PARENT_SCOPE)
