# Source for special functions
set(SRC_SPECIAL
    ${CMAKE_CURRENT_SOURCE_DIR}/spherical_harmonics.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/circular_harmonics.f90)

set(SRC ${SRC} ${SRC_SPECIAL} PARENT_SCOPE)