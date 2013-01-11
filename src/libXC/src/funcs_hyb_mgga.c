#include "util.h"

extern XC(func_info_type) XC(func_info_hyb_mgga_xc_m05);
extern XC(func_info_type) XC(func_info_hyb_mgga_xc_m05_2x);
extern XC(func_info_type) XC(func_info_hyb_mgga_xc_b88b95);
extern XC(func_info_type) XC(func_info_hyb_mgga_xc_b86b95);
extern XC(func_info_type) XC(func_info_hyb_mgga_xc_pw86b95);
extern XC(func_info_type) XC(func_info_hyb_mgga_xc_bb1k);
extern XC(func_info_type) XC(func_info_hyb_mgga_xc_m06_hf);
extern XC(func_info_type) XC(func_info_hyb_mgga_xc_mpw1b95);
extern XC(func_info_type) XC(func_info_hyb_mgga_xc_mpwb1k);
extern XC(func_info_type) XC(func_info_hyb_mgga_xc_x1b95);
extern XC(func_info_type) XC(func_info_hyb_mgga_xc_xb1k);
extern XC(func_info_type) XC(func_info_hyb_mgga_xc_m06);
extern XC(func_info_type) XC(func_info_hyb_mgga_xc_m06_2x);


const XC(func_info_type) *XC(hyb_mgga_known_funct)[] = {
  &XC(func_info_hyb_mgga_xc_m05),
  &XC(func_info_hyb_mgga_xc_m05_2x),
  &XC(func_info_hyb_mgga_xc_b88b95),
  &XC(func_info_hyb_mgga_xc_b86b95),
  &XC(func_info_hyb_mgga_xc_pw86b95),
  &XC(func_info_hyb_mgga_xc_bb1k),
  &XC(func_info_hyb_mgga_xc_m06_hf),
  &XC(func_info_hyb_mgga_xc_mpw1b95),
  &XC(func_info_hyb_mgga_xc_mpwb1k),
  &XC(func_info_hyb_mgga_xc_x1b95),
  &XC(func_info_hyb_mgga_xc_xb1k),
  &XC(func_info_hyb_mgga_xc_m06),
  &XC(func_info_hyb_mgga_xc_m06_2x),
  NULL
};
