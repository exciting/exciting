# Library for the cell
set(SRC_CRYSTAL
    ${CMAKE_CURRENT_SOURCE_DIR}/crystal_cell.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/crystal_symmetry.f90)

set(SRC ${SRC} ${SRC_CRYSTAL} PARENT_SCOPE)