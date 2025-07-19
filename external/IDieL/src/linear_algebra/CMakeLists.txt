# Computation of auxiliary vectors
set(SRC_LINALG ${SRC_LINALG}
    ${CMAKE_CURRENT_SOURCE_DIR}/linalg.f90)

set(SRC ${SRC} ${SRC_LINALG} PARENT_SCOPE)
