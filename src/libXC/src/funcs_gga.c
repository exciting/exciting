#include "util.h"

extern XC(func_info_type) XC(func_info_gga_x_pbe);
extern XC(func_info_type) XC(func_info_gga_x_pbe_r);
extern XC(func_info_type) XC(func_info_gga_x_b86);
extern XC(func_info_type) XC(func_info_gga_x_b86_r);
extern XC(func_info_type) XC(func_info_gga_x_b86_mgc);
extern XC(func_info_type) XC(func_info_gga_x_b88);
extern XC(func_info_type) XC(func_info_gga_x_g96);
extern XC(func_info_type) XC(func_info_gga_x_pw86);
extern XC(func_info_type) XC(func_info_gga_x_pw91);
extern XC(func_info_type) XC(func_info_gga_x_optx);
extern XC(func_info_type) XC(func_info_gga_x_dk87_r1);
extern XC(func_info_type) XC(func_info_gga_x_dk87_r2);
extern XC(func_info_type) XC(func_info_gga_x_lg93);
extern XC(func_info_type) XC(func_info_gga_x_ft97_a);
extern XC(func_info_type) XC(func_info_gga_x_ft97_b);
extern XC(func_info_type) XC(func_info_gga_x_pbe_sol);
extern XC(func_info_type) XC(func_info_gga_x_rpbe);
extern XC(func_info_type) XC(func_info_gga_x_wc);
extern XC(func_info_type) XC(func_info_gga_x_mpw91);
extern XC(func_info_type) XC(func_info_gga_x_am05);
extern XC(func_info_type) XC(func_info_gga_x_pbea);
extern XC(func_info_type) XC(func_info_gga_x_mpbe);
extern XC(func_info_type) XC(func_info_gga_x_xpbe);
extern XC(func_info_type) XC(func_info_gga_x_2d_b86_mgc);
extern XC(func_info_type) XC(func_info_gga_x_bayesian);
extern XC(func_info_type) XC(func_info_gga_x_pbe_jsjr);
extern XC(func_info_type) XC(func_info_gga_x_2d_b88);
extern XC(func_info_type) XC(func_info_gga_x_2d_b86);
extern XC(func_info_type) XC(func_info_gga_x_2d_pbe);
extern XC(func_info_type) XC(func_info_gga_c_pbe);
extern XC(func_info_type) XC(func_info_gga_c_lyp);
extern XC(func_info_type) XC(func_info_gga_c_p86);
extern XC(func_info_type) XC(func_info_gga_c_pbe_sol);
extern XC(func_info_type) XC(func_info_gga_c_pw91);
extern XC(func_info_type) XC(func_info_gga_c_am05);
extern XC(func_info_type) XC(func_info_gga_c_xpbe);
extern XC(func_info_type) XC(func_info_gga_c_lm);
extern XC(func_info_type) XC(func_info_gga_c_pbe_jrgx);
extern XC(func_info_type) XC(func_info_gga_x_optb88_vdw);
extern XC(func_info_type) XC(func_info_gga_x_pbek1_vdw);
extern XC(func_info_type) XC(func_info_gga_x_optpbe_vdw);
extern XC(func_info_type) XC(func_info_gga_x_rge2);
extern XC(func_info_type) XC(func_info_gga_c_rge2);
extern XC(func_info_type) XC(func_info_gga_xc_lb);
extern XC(func_info_type) XC(func_info_gga_xc_hcth_93);
extern XC(func_info_type) XC(func_info_gga_xc_hcth_120);
extern XC(func_info_type) XC(func_info_gga_xc_hcth_147);
extern XC(func_info_type) XC(func_info_gga_xc_hcth_407);
extern XC(func_info_type) XC(func_info_gga_xc_edf1);
extern XC(func_info_type) XC(func_info_gga_xc_xlyp);
extern XC(func_info_type) XC(func_info_gga_xc_b97);
extern XC(func_info_type) XC(func_info_gga_xc_b97_1);
extern XC(func_info_type) XC(func_info_gga_xc_b97_2);
extern XC(func_info_type) XC(func_info_gga_xc_b97_d);
extern XC(func_info_type) XC(func_info_gga_xc_b97_k);
extern XC(func_info_type) XC(func_info_gga_xc_b97_3);
extern XC(func_info_type) XC(func_info_gga_xc_pbe1w);
extern XC(func_info_type) XC(func_info_gga_xc_mpwlyp1w);
extern XC(func_info_type) XC(func_info_gga_xc_pbelyp1w);
extern XC(func_info_type) XC(func_info_gga_xc_sb98_1a);
extern XC(func_info_type) XC(func_info_gga_xc_sb98_1b);
extern XC(func_info_type) XC(func_info_gga_xc_sb98_1c);
extern XC(func_info_type) XC(func_info_gga_xc_sb98_2a);
extern XC(func_info_type) XC(func_info_gga_xc_sb98_2b);
extern XC(func_info_type) XC(func_info_gga_xc_sb98_2c);


const XC(func_info_type) *XC(gga_known_funct)[] = {
  &XC(func_info_gga_x_pbe),
  &XC(func_info_gga_x_pbe_r),
  &XC(func_info_gga_x_b86),
  &XC(func_info_gga_x_b86_r),
  &XC(func_info_gga_x_b86_mgc),
  &XC(func_info_gga_x_b88),
  &XC(func_info_gga_x_g96),
  &XC(func_info_gga_x_pw86),
  &XC(func_info_gga_x_pw91),
  &XC(func_info_gga_x_optx),
  &XC(func_info_gga_x_dk87_r1),
  &XC(func_info_gga_x_dk87_r2),
  &XC(func_info_gga_x_lg93),
  &XC(func_info_gga_x_ft97_a),
  &XC(func_info_gga_x_ft97_b),
  &XC(func_info_gga_x_pbe_sol),
  &XC(func_info_gga_x_rpbe),
  &XC(func_info_gga_x_wc),
  &XC(func_info_gga_x_mpw91),
  &XC(func_info_gga_x_am05),
  &XC(func_info_gga_x_pbea),
  &XC(func_info_gga_x_mpbe),
  &XC(func_info_gga_x_xpbe),
  &XC(func_info_gga_x_2d_b86_mgc),
  &XC(func_info_gga_x_bayesian),
  &XC(func_info_gga_x_pbe_jsjr),
  &XC(func_info_gga_x_2d_b88),
  &XC(func_info_gga_x_2d_b86),
  &XC(func_info_gga_x_2d_pbe),
  &XC(func_info_gga_c_pbe),
  &XC(func_info_gga_c_lyp),
  &XC(func_info_gga_c_p86),
  &XC(func_info_gga_c_pbe_sol),
  &XC(func_info_gga_c_pw91),
  &XC(func_info_gga_c_am05),
  &XC(func_info_gga_c_xpbe),
  &XC(func_info_gga_c_lm),
  &XC(func_info_gga_c_pbe_jrgx),
  &XC(func_info_gga_x_optb88_vdw),
  &XC(func_info_gga_x_pbek1_vdw),
  &XC(func_info_gga_x_optpbe_vdw),
  &XC(func_info_gga_x_rge2),
  &XC(func_info_gga_c_rge2),
  &XC(func_info_gga_xc_lb),
  &XC(func_info_gga_xc_hcth_93),
  &XC(func_info_gga_xc_hcth_120),
  &XC(func_info_gga_xc_hcth_147),
  &XC(func_info_gga_xc_hcth_407),
  &XC(func_info_gga_xc_edf1),
  &XC(func_info_gga_xc_xlyp),
  &XC(func_info_gga_xc_b97),
  &XC(func_info_gga_xc_b97_1),
  &XC(func_info_gga_xc_b97_2),
  &XC(func_info_gga_xc_b97_d),
  &XC(func_info_gga_xc_b97_k),
  &XC(func_info_gga_xc_b97_3),
  &XC(func_info_gga_xc_pbe1w),
  &XC(func_info_gga_xc_mpwlyp1w),
  &XC(func_info_gga_xc_pbelyp1w),
  &XC(func_info_gga_xc_sb98_1a),
  &XC(func_info_gga_xc_sb98_1b),
  &XC(func_info_gga_xc_sb98_1c),
  &XC(func_info_gga_xc_sb98_2a),
  &XC(func_info_gga_xc_sb98_2b),
  &XC(func_info_gga_xc_sb98_2c),
  NULL
};
