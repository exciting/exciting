#include "util.h"

extern xc_func_info_type xc_func_info_hyb_lda_xc_lda0;
extern xc_func_info_type xc_func_info_hyb_lda_xc_cam_lda0;
extern xc_func_info_type xc_func_info_hyb_lda_xc_bn05;
extern xc_func_info_type xc_func_info_hyb_lda_x_erf;

const xc_func_info_type *xc_hyb_lda_known_funct[] = {
  &xc_func_info_hyb_lda_xc_lda0,
  &xc_func_info_hyb_lda_xc_cam_lda0,
  &xc_func_info_hyb_lda_xc_bn05,
  &xc_func_info_hyb_lda_x_erf,
  NULL
};
