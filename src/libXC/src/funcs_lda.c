#include "util.h"

extern XC(func_info_type) XC(func_info_lda_x);
extern XC(func_info_type) XC(func_info_lda_c_wigner);
extern XC(func_info_type) XC(func_info_lda_c_rpa);
extern XC(func_info_type) XC(func_info_lda_c_hl);
extern XC(func_info_type) XC(func_info_lda_c_gl);
extern XC(func_info_type) XC(func_info_lda_c_xalpha);
extern XC(func_info_type) XC(func_info_lda_c_vwn);
extern XC(func_info_type) XC(func_info_lda_c_vwn_rpa);
extern XC(func_info_type) XC(func_info_lda_c_pz);
extern XC(func_info_type) XC(func_info_lda_c_pz_mod);
extern XC(func_info_type) XC(func_info_lda_c_ob_pz);
extern XC(func_info_type) XC(func_info_lda_c_pw);
extern XC(func_info_type) XC(func_info_lda_c_pw_mod);
extern XC(func_info_type) XC(func_info_lda_c_ob_pw);
extern XC(func_info_type) XC(func_info_lda_c_2d_amgb);
extern XC(func_info_type) XC(func_info_lda_c_2d_prm);
extern XC(func_info_type) XC(func_info_lda_c_vbh);
extern XC(func_info_type) XC(func_info_lda_c_1d_csc);
extern XC(func_info_type) XC(func_info_lda_x_2d);
extern XC(func_info_type) XC(func_info_lda_xc_teter93);
extern XC(func_info_type) XC(func_info_lda_x_1d);
extern XC(func_info_type) XC(func_info_lda_c_ml1);
extern XC(func_info_type) XC(func_info_lda_c_ml2);
extern XC(func_info_type) XC(func_info_lda_c_gombas);
extern XC(func_info_type) XC(func_info_lda_c_pw_rpa);
extern XC(func_info_type) XC(func_info_lda_c_1d_loos);
extern XC(func_info_type) XC(func_info_lda_c_rc04);
extern XC(func_info_type) XC(func_info_lda_c_vwn_1);
extern XC(func_info_type) XC(func_info_lda_c_vwn_2);
extern XC(func_info_type) XC(func_info_lda_c_vwn_3);
extern XC(func_info_type) XC(func_info_lda_c_vwn_4);
extern XC(func_info_type) XC(func_info_lda_k_tf);
extern XC(func_info_type) XC(func_info_lda_k_lp);


const XC(func_info_type) *XC(lda_known_funct)[] = {
  &XC(func_info_lda_x),
  &XC(func_info_lda_c_wigner),
  &XC(func_info_lda_c_rpa),
  &XC(func_info_lda_c_hl),
  &XC(func_info_lda_c_gl),
  &XC(func_info_lda_c_xalpha),
  &XC(func_info_lda_c_vwn),
  &XC(func_info_lda_c_vwn_rpa),
  &XC(func_info_lda_c_pz),
  &XC(func_info_lda_c_pz_mod),
  &XC(func_info_lda_c_ob_pz),
  &XC(func_info_lda_c_pw),
  &XC(func_info_lda_c_pw_mod),
  &XC(func_info_lda_c_ob_pw),
  &XC(func_info_lda_c_2d_amgb),
  &XC(func_info_lda_c_2d_prm),
  &XC(func_info_lda_c_vbh),
  &XC(func_info_lda_c_1d_csc),
  &XC(func_info_lda_x_2d),
  &XC(func_info_lda_xc_teter93),
  &XC(func_info_lda_x_1d),
  &XC(func_info_lda_c_ml1),
  &XC(func_info_lda_c_ml2),
  &XC(func_info_lda_c_gombas),
  &XC(func_info_lda_c_pw_rpa),
  &XC(func_info_lda_c_1d_loos),
  &XC(func_info_lda_c_rc04),
  &XC(func_info_lda_c_vwn_1),
  &XC(func_info_lda_c_vwn_2),
  &XC(func_info_lda_c_vwn_3),
  &XC(func_info_lda_c_vwn_4),
  &XC(func_info_lda_k_tf),
  &XC(func_info_lda_k_lp),
  NULL
};
