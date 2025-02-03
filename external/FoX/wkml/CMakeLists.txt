# m_wkml_coverage.f90: m_wkml_coverage.m4
#         m4 -I ../m4 $< > $@
IF(NOT ${M4_PROGRAM} MATCHES "NOTFOUND")
  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/m_wkml_coverage.F90
    COMMAND ${M4_PROGRAM} -I ../m4 ${CMAKE_CURRENT_SOURCE_DIR}/m_wkml_coverage.m4 > ${CMAKE_CURRENT_SOURCE_DIR}/m_wkml_coverage.F90
    DEPENDS m_wkml_coverage.m4
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
ENDIF(NOT ${M4_PROGRAM} MATCHES "NOTFOUND")

set(fox_wkml_srcs
  FoX_wkml.f90
)

set(fox_wkml_srcs_pp
  m_wkml_lowlevel.F90
  m_wkml_color.F90
  m_wkml_styling.F90
  m_wkml_features.F90
  m_wkml_coverage.F90
  m_wkml_core.F90
  m_wkml_contours.F90
  m_contours.F90
  m_wkml_color_def.F90
  m_wkml_chart.F90
)

MK_LIBS(wkml utils wxml)