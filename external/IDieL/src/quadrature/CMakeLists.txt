# Library for the spherical quadrature
set(SRC_QUADRATURE 
    ${CMAKE_CURRENT_SOURCE_DIR}/sph_quadrature.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/circle_quadrature.f90)

set(SRC ${SRC} ${SRC_QUADRATURE} PARENT_SCOPE)