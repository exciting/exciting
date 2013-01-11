#include "util.h"

extern XC(func_info_type) XC(func_info_hyb_gga_xc_b3pw91);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_b3lyp);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_b3p86);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_o3lyp);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_mpw1k);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_pbeh);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_b97);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_b97_1);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_b97_2);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_x3lyp);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_b1wc);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_b97_k);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_b97_3);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_mpw3pw);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_b1lyp);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_b1pw91);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_mpw1pw);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_mpw3lyp);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_sb98_1a);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_sb98_1b);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_sb98_1c);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_sb98_2a);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_sb98_2b);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_sb98_2c);
extern XC(func_info_type) XC(func_info_hyb_gga_x_sogga11_x);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_hse03);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_hse06);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_hjs_pbe);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_hjs_pbe_sol);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_hjs_b88);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_hjs_b97x);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_cam_b3lyp);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_tuned_cam_b3lyp);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_bhandh);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_bhandhlyp);
extern XC(func_info_type) XC(func_info_hyb_gga_xc_mb3lyp_rc04);


const XC(func_info_type) *XC(hyb_gga_known_funct)[] = {
  &XC(func_info_hyb_gga_xc_b3pw91),
  &XC(func_info_hyb_gga_xc_b3lyp),
  &XC(func_info_hyb_gga_xc_b3p86),
  &XC(func_info_hyb_gga_xc_o3lyp),
  &XC(func_info_hyb_gga_xc_mpw1k),
  &XC(func_info_hyb_gga_xc_pbeh),
  &XC(func_info_hyb_gga_xc_b97),
  &XC(func_info_hyb_gga_xc_b97_1),
  &XC(func_info_hyb_gga_xc_b97_2),
  &XC(func_info_hyb_gga_xc_x3lyp),
  &XC(func_info_hyb_gga_xc_b1wc),
  &XC(func_info_hyb_gga_xc_b97_k),
  &XC(func_info_hyb_gga_xc_b97_3),
  &XC(func_info_hyb_gga_xc_mpw3pw),
  &XC(func_info_hyb_gga_xc_b1lyp),
  &XC(func_info_hyb_gga_xc_b1pw91),
  &XC(func_info_hyb_gga_xc_mpw1pw),
  &XC(func_info_hyb_gga_xc_mpw3lyp),
  &XC(func_info_hyb_gga_xc_sb98_1a),
  &XC(func_info_hyb_gga_xc_sb98_1b),
  &XC(func_info_hyb_gga_xc_sb98_1c),
  &XC(func_info_hyb_gga_xc_sb98_2a),
  &XC(func_info_hyb_gga_xc_sb98_2b),
  &XC(func_info_hyb_gga_xc_sb98_2c),
  &XC(func_info_hyb_gga_x_sogga11_x),
  &XC(func_info_hyb_gga_xc_hse03),
  &XC(func_info_hyb_gga_xc_hse06),
  &XC(func_info_hyb_gga_xc_hjs_pbe),
  &XC(func_info_hyb_gga_xc_hjs_pbe_sol),
  &XC(func_info_hyb_gga_xc_hjs_b88),
  &XC(func_info_hyb_gga_xc_hjs_b97x),
  &XC(func_info_hyb_gga_xc_cam_b3lyp),
  &XC(func_info_hyb_gga_xc_tuned_cam_b3lyp),
  &XC(func_info_hyb_gga_xc_bhandh),
  &XC(func_info_hyb_gga_xc_bhandhlyp),
  &XC(func_info_hyb_gga_xc_mb3lyp_rc04),
  NULL
};
