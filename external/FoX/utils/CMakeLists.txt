set(fox_utils_srcs
  FoX_utils.f90
)

set(fox_utils_srcs_pp
  fox_m_utils_mtprng.F90
  fox_m_utils_uri.F90
  fox_m_utils_uuid.F90
)

MK_LIBS(utils fsys)