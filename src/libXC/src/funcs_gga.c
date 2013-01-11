#include "util.h"

extern XC(func_info_type) XC(func_info_gga_c_op_xalpha);
extern XC(func_info_type) XC(func_info_gga_c_op_g96);
extern XC(func_info_type) XC(func_info_gga_c_op_pbe);
extern XC(func_info_type) XC(func_info_gga_c_op_b88);
extern XC(func_info_type) XC(func_info_gga_c_ft97);
extern XC(func_info_type) XC(func_info_gga_c_spbe);
extern XC(func_info_type) XC(func_info_gga_x_ssb_sw);
extern XC(func_info_type) XC(func_info_gga_x_ssb);
extern XC(func_info_type) XC(func_info_gga_x_ssb_d);
extern XC(func_info_type) XC(func_info_gga_xc_hcth_407p);
extern XC(func_info_type) XC(func_info_gga_xc_hcth_p76);
extern XC(func_info_type) XC(func_info_gga_xc_hcth_p14);
extern XC(func_info_type) XC(func_info_gga_xc_b97_gga1);
extern XC(func_info_type) XC(func_info_gga_xc_hcth_a);
extern XC(func_info_type) XC(func_info_gga_x_bpccac);
extern XC(func_info_type) XC(func_info_gga_c_revtca);
extern XC(func_info_type) XC(func_info_gga_c_tca);
extern XC(func_info_type) XC(func_info_gga_x_pbe);
extern XC(func_info_type) XC(func_info_gga_x_pbe_r);
extern XC(func_info_type) XC(func_info_gga_x_b86);
extern XC(func_info_type) XC(func_info_gga_x_herman);
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
extern XC(func_info_type) XC(func_info_gga_x_rpw86);
extern XC(func_info_type) XC(func_info_gga_x_kt1);
extern XC(func_info_type) XC(func_info_gga_xc_kt2);
extern XC(func_info_type) XC(func_info_gga_c_wl);
extern XC(func_info_type) XC(func_info_gga_c_wi);
extern XC(func_info_type) XC(func_info_gga_x_mb88);
extern XC(func_info_type) XC(func_info_gga_x_sogga);
extern XC(func_info_type) XC(func_info_gga_x_sogga11);
extern XC(func_info_type) XC(func_info_gga_c_sogga11);
extern XC(func_info_type) XC(func_info_gga_c_wi0);
extern XC(func_info_type) XC(func_info_gga_xc_th1);
extern XC(func_info_type) XC(func_info_gga_xc_th2);
extern XC(func_info_type) XC(func_info_gga_xc_th3);
extern XC(func_info_type) XC(func_info_gga_xc_th4);
extern XC(func_info_type) XC(func_info_gga_x_c09x);
extern XC(func_info_type) XC(func_info_gga_c_sogga11_x);
extern XC(func_info_type) XC(func_info_gga_x_lb);
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
extern XC(func_info_type) XC(func_info_gga_x_lbm);
extern XC(func_info_type) XC(func_info_gga_x_ol2);
extern XC(func_info_type) XC(func_info_gga_x_apbe);
extern XC(func_info_type) XC(func_info_gga_k_apbe);
extern XC(func_info_type) XC(func_info_gga_c_apbe);
extern XC(func_info_type) XC(func_info_gga_k_tw1);
extern XC(func_info_type) XC(func_info_gga_k_tw2);
extern XC(func_info_type) XC(func_info_gga_k_tw3);
extern XC(func_info_type) XC(func_info_gga_k_tw4);
extern XC(func_info_type) XC(func_info_gga_x_htbs);
extern XC(func_info_type) XC(func_info_gga_x_airy);
extern XC(func_info_type) XC(func_info_gga_x_lag);
extern XC(func_info_type) XC(func_info_gga_xc_mohlyp);
extern XC(func_info_type) XC(func_info_gga_xc_mohlyp2);
extern XC(func_info_type) XC(func_info_gga_xc_th_fl);
extern XC(func_info_type) XC(func_info_gga_xc_th_fc);
extern XC(func_info_type) XC(func_info_gga_xc_th_fcfo);
extern XC(func_info_type) XC(func_info_gga_xc_th_fco);
extern XC(func_info_type) XC(func_info_gga_c_optc);
extern XC(func_info_type) XC(func_info_gga_k_vw);
extern XC(func_info_type) XC(func_info_gga_k_ge2);
extern XC(func_info_type) XC(func_info_gga_k_golden);
extern XC(func_info_type) XC(func_info_gga_k_yt65);
extern XC(func_info_type) XC(func_info_gga_k_baltin);
extern XC(func_info_type) XC(func_info_gga_k_lieb);
extern XC(func_info_type) XC(func_info_gga_k_absr1);
extern XC(func_info_type) XC(func_info_gga_k_absr2);
extern XC(func_info_type) XC(func_info_gga_k_gr);
extern XC(func_info_type) XC(func_info_gga_k_ludena);
extern XC(func_info_type) XC(func_info_gga_k_gp85);
extern XC(func_info_type) XC(func_info_gga_k_pearson);
extern XC(func_info_type) XC(func_info_gga_k_ol1);
extern XC(func_info_type) XC(func_info_gga_k_ol2);
extern XC(func_info_type) XC(func_info_gga_k_fr_b88);
extern XC(func_info_type) XC(func_info_gga_k_fr_pw86);
extern XC(func_info_type) XC(func_info_gga_k_dk);
extern XC(func_info_type) XC(func_info_gga_k_perdew);
extern XC(func_info_type) XC(func_info_gga_k_vsk);
extern XC(func_info_type) XC(func_info_gga_k_vjks);
extern XC(func_info_type) XC(func_info_gga_k_ernzerhof);
extern XC(func_info_type) XC(func_info_gga_k_lc94);
extern XC(func_info_type) XC(func_info_gga_k_llp);
extern XC(func_info_type) XC(func_info_gga_k_thakkar);
extern XC(func_info_type) XC(func_info_gga_x_wpbeh);
extern XC(func_info_type) XC(func_info_gga_x_hjs_pbe);
extern XC(func_info_type) XC(func_info_gga_x_hjs_pbe_sol);
extern XC(func_info_type) XC(func_info_gga_x_hjs_b88);
extern XC(func_info_type) XC(func_info_gga_x_hjs_b97x);
extern XC(func_info_type) XC(func_info_gga_x_ityh);


const XC(func_info_type) *XC(gga_known_funct)[] = {
  &XC(func_info_gga_c_op_xalpha),
  &XC(func_info_gga_c_op_g96),
  &XC(func_info_gga_c_op_pbe),
  &XC(func_info_gga_c_op_b88),
  &XC(func_info_gga_c_ft97),
  &XC(func_info_gga_c_spbe),
  &XC(func_info_gga_x_ssb_sw),
  &XC(func_info_gga_x_ssb),
  &XC(func_info_gga_x_ssb_d),
  &XC(func_info_gga_xc_hcth_407p),
  &XC(func_info_gga_xc_hcth_p76),
  &XC(func_info_gga_xc_hcth_p14),
  &XC(func_info_gga_xc_b97_gga1),
  &XC(func_info_gga_xc_hcth_a),
  &XC(func_info_gga_x_bpccac),
  &XC(func_info_gga_c_revtca),
  &XC(func_info_gga_c_tca),
  &XC(func_info_gga_x_pbe),
  &XC(func_info_gga_x_pbe_r),
  &XC(func_info_gga_x_b86),
  &XC(func_info_gga_x_herman),
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
  &XC(func_info_gga_x_rpw86),
  &XC(func_info_gga_x_kt1),
  &XC(func_info_gga_xc_kt2),
  &XC(func_info_gga_c_wl),
  &XC(func_info_gga_c_wi),
  &XC(func_info_gga_x_mb88),
  &XC(func_info_gga_x_sogga),
  &XC(func_info_gga_x_sogga11),
  &XC(func_info_gga_c_sogga11),
  &XC(func_info_gga_c_wi0),
  &XC(func_info_gga_xc_th1),
  &XC(func_info_gga_xc_th2),
  &XC(func_info_gga_xc_th3),
  &XC(func_info_gga_xc_th4),
  &XC(func_info_gga_x_c09x),
  &XC(func_info_gga_c_sogga11_x),
  &XC(func_info_gga_x_lb),
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
  &XC(func_info_gga_x_lbm),
  &XC(func_info_gga_x_ol2),
  &XC(func_info_gga_x_apbe),
  &XC(func_info_gga_k_apbe),
  &XC(func_info_gga_c_apbe),
  &XC(func_info_gga_k_tw1),
  &XC(func_info_gga_k_tw2),
  &XC(func_info_gga_k_tw3),
  &XC(func_info_gga_k_tw4),
  &XC(func_info_gga_x_htbs),
  &XC(func_info_gga_x_airy),
  &XC(func_info_gga_x_lag),
  &XC(func_info_gga_xc_mohlyp),
  &XC(func_info_gga_xc_mohlyp2),
  &XC(func_info_gga_xc_th_fl),
  &XC(func_info_gga_xc_th_fc),
  &XC(func_info_gga_xc_th_fcfo),
  &XC(func_info_gga_xc_th_fco),
  &XC(func_info_gga_c_optc),
  &XC(func_info_gga_k_vw),
  &XC(func_info_gga_k_ge2),
  &XC(func_info_gga_k_golden),
  &XC(func_info_gga_k_yt65),
  &XC(func_info_gga_k_baltin),
  &XC(func_info_gga_k_lieb),
  &XC(func_info_gga_k_absr1),
  &XC(func_info_gga_k_absr2),
  &XC(func_info_gga_k_gr),
  &XC(func_info_gga_k_ludena),
  &XC(func_info_gga_k_gp85),
  &XC(func_info_gga_k_pearson),
  &XC(func_info_gga_k_ol1),
  &XC(func_info_gga_k_ol2),
  &XC(func_info_gga_k_fr_b88),
  &XC(func_info_gga_k_fr_pw86),
  &XC(func_info_gga_k_dk),
  &XC(func_info_gga_k_perdew),
  &XC(func_info_gga_k_vsk),
  &XC(func_info_gga_k_vjks),
  &XC(func_info_gga_k_ernzerhof),
  &XC(func_info_gga_k_lc94),
  &XC(func_info_gga_k_llp),
  &XC(func_info_gga_k_thakkar),
  &XC(func_info_gga_x_wpbeh),
  &XC(func_info_gga_x_hjs_pbe),
  &XC(func_info_gga_x_hjs_pbe_sol),
  &XC(func_info_gga_x_hjs_b88),
  &XC(func_info_gga_x_hjs_b97x),
  &XC(func_info_gga_x_ityh),
  NULL
};
