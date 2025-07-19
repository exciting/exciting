# Dielectric averaging functions
set(SRC_DIELECTRIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/idiel.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/definitions/idiel_common.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/definitions/idiel_2d.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/definitions/idiel_3d.f90)

set(SRC ${SRC} ${SRC_DIELECTRIC} PARENT_SCOPE)