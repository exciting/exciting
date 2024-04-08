# Local components
if(NOT CPUBACKEND)
    set(SRC_LINALG ${CMAKE_CURRENT_SOURCE_DIR}/device/linalg_device_common.f90)
else()
    set(SRC_LINALG ${CMAKE_CURRENT_SOURCE_DIR}/host/linalg_device_common.f90)
endif()
# Push back to SRC
set(SRC ${SRC} ${SRC_LINALG} PARENT_SCOPE)
