# Source for special functions
set(SRC_SPLINE
    ${CMAKE_CURRENT_SOURCE_DIR}/cubic_spline_1d.f90)

set(SRC ${SRC} ${SRC_SPLINE} PARENT_SCOPE)