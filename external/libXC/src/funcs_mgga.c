#include "util.h"

extern xc_func_info_type xc_func_info_mgga_c_dldf;
extern xc_func_info_type xc_func_info_mgga_xc_zlp;
extern xc_func_info_type xc_func_info_mgga_xc_otpss_d;
extern xc_func_info_type xc_func_info_mgga_c_cs;
extern xc_func_info_type xc_func_info_mgga_c_mn12_sx;
extern xc_func_info_type xc_func_info_mgga_c_mn12_l;
extern xc_func_info_type xc_func_info_mgga_c_m11_l;
extern xc_func_info_type xc_func_info_mgga_c_m11;
extern xc_func_info_type xc_func_info_mgga_c_m08_so;
extern xc_func_info_type xc_func_info_mgga_c_m08_hx;
extern xc_func_info_type xc_func_info_mgga_c_revm11;
extern xc_func_info_type xc_func_info_mgga_x_lta;
extern xc_func_info_type xc_func_info_mgga_x_tpss;
extern xc_func_info_type xc_func_info_mgga_x_m06_l;
extern xc_func_info_type xc_func_info_mgga_x_gvt4;
extern xc_func_info_type xc_func_info_mgga_x_tau_hcth;
extern xc_func_info_type xc_func_info_mgga_x_br89;
extern xc_func_info_type xc_func_info_mgga_x_bj06;
extern xc_func_info_type xc_func_info_mgga_x_tb09;
extern xc_func_info_type xc_func_info_mgga_x_rpp09;
extern xc_func_info_type xc_func_info_mgga_x_2d_prhg07;
extern xc_func_info_type xc_func_info_mgga_x_2d_prhg07_prp10;
extern xc_func_info_type xc_func_info_mgga_x_revtpss;
extern xc_func_info_type xc_func_info_mgga_x_pkzb;
extern xc_func_info_type xc_func_info_mgga_x_br89_1;
extern xc_func_info_type xc_func_info_mgga_k_pgsl025;
extern xc_func_info_type xc_func_info_mgga_x_ms0;
extern xc_func_info_type xc_func_info_mgga_x_ms1;
extern xc_func_info_type xc_func_info_mgga_x_ms2;
extern xc_func_info_type xc_func_info_mgga_x_th;
extern xc_func_info_type xc_func_info_mgga_x_m11_l;
extern xc_func_info_type xc_func_info_mgga_x_mn12_l;
extern xc_func_info_type xc_func_info_mgga_x_ms2_rev;
extern xc_func_info_type xc_func_info_mgga_xc_cc06;
extern xc_func_info_type xc_func_info_mgga_x_mk00;
extern xc_func_info_type xc_func_info_mgga_c_tpss;
extern xc_func_info_type xc_func_info_mgga_c_vsxc;
extern xc_func_info_type xc_func_info_mgga_c_m06_l;
extern xc_func_info_type xc_func_info_mgga_c_m06_hf;
extern xc_func_info_type xc_func_info_mgga_c_m06;
extern xc_func_info_type xc_func_info_mgga_c_m06_2x;
extern xc_func_info_type xc_func_info_mgga_c_m05;
extern xc_func_info_type xc_func_info_mgga_c_m05_2x;
extern xc_func_info_type xc_func_info_mgga_c_pkzb;
extern xc_func_info_type xc_func_info_mgga_c_bc95;
extern xc_func_info_type xc_func_info_mgga_c_revtpss;
extern xc_func_info_type xc_func_info_mgga_xc_tpsslyp1w;
extern xc_func_info_type xc_func_info_mgga_x_mk00b;
extern xc_func_info_type xc_func_info_mgga_x_bloc;
extern xc_func_info_type xc_func_info_mgga_x_modtpss;
extern xc_func_info_type xc_func_info_mgga_c_tpssloc;
extern xc_func_info_type xc_func_info_mgga_x_mbeef;
extern xc_func_info_type xc_func_info_mgga_x_mbeefvdw;
extern xc_func_info_type xc_func_info_mgga_c_tm;
extern xc_func_info_type xc_func_info_mgga_xc_b97m_v;
extern xc_func_info_type xc_func_info_mgga_x_jk;
extern xc_func_info_type xc_func_info_mgga_x_mvs;
extern xc_func_info_type xc_func_info_mgga_x_mn15_l;
extern xc_func_info_type xc_func_info_mgga_c_mn15_l;
extern xc_func_info_type xc_func_info_mgga_x_scan;
extern xc_func_info_type xc_func_info_mgga_c_scan;
extern xc_func_info_type xc_func_info_mgga_c_mn15;
extern xc_func_info_type xc_func_info_mgga_x_b00;
extern xc_func_info_type xc_func_info_mgga_xc_hle17;
extern xc_func_info_type xc_func_info_mgga_c_scan_rvv10;
extern xc_func_info_type xc_func_info_mgga_x_revm06_l;
extern xc_func_info_type xc_func_info_mgga_c_revm06_l;
extern xc_func_info_type xc_func_info_mgga_x_rtpss;
extern xc_func_info_type xc_func_info_mgga_x_ms2b;
extern xc_func_info_type xc_func_info_mgga_x_ms2bs;
extern xc_func_info_type xc_func_info_mgga_x_mvsb;
extern xc_func_info_type xc_func_info_mgga_x_mvsbs;
extern xc_func_info_type xc_func_info_mgga_c_revm06;
extern xc_func_info_type xc_func_info_mgga_c_m06_sx;
extern xc_func_info_type xc_func_info_mgga_x_ft98;
extern xc_func_info_type xc_func_info_mgga_c_tpss_gaussian;
extern xc_func_info_type xc_func_info_mgga_x_eel;
extern xc_func_info_type xc_func_info_mgga_c_cf22d;
extern xc_func_info_type xc_func_info_mgga_x_lak;
extern xc_func_info_type xc_func_info_mgga_c_cc;
extern xc_func_info_type xc_func_info_mgga_c_ccalda;
extern xc_func_info_type xc_func_info_mgga_c_rregtm;
extern xc_func_info_type xc_func_info_mgga_c_b94;
extern xc_func_info_type xc_func_info_mgga_x_rscan;
extern xc_func_info_type xc_func_info_mgga_c_rscan;
extern xc_func_info_type xc_func_info_mgga_x_r2scan;
extern xc_func_info_type xc_func_info_mgga_c_r2scan;
extern xc_func_info_type xc_func_info_mgga_x_tm;
extern xc_func_info_type xc_func_info_mgga_x_vt84;
extern xc_func_info_type xc_func_info_mgga_x_sa_tpss;
extern xc_func_info_type xc_func_info_mgga_k_pc07;
extern xc_func_info_type xc_func_info_mgga_c_kcis;
extern xc_func_info_type xc_func_info_mgga_xc_lp90;
extern xc_func_info_type xc_func_info_mgga_c_b88;
extern xc_func_info_type xc_func_info_mgga_x_gx;
extern xc_func_info_type xc_func_info_mgga_x_pbe_gx;
extern xc_func_info_type xc_func_info_mgga_x_revscan;
extern xc_func_info_type xc_func_info_mgga_c_revscan;
extern xc_func_info_type xc_func_info_mgga_c_scan_vv10;
extern xc_func_info_type xc_func_info_mgga_c_revscan_vv10;
extern xc_func_info_type xc_func_info_mgga_x_br89_explicit;
extern xc_func_info_type xc_func_info_mgga_x_br89_explicit_1;
extern xc_func_info_type xc_func_info_mgga_x_regtpss;
extern xc_func_info_type xc_func_info_mgga_x_2d_js17;
extern xc_func_info_type xc_func_info_mgga_k_l04;
extern xc_func_info_type xc_func_info_mgga_k_l06;
extern xc_func_info_type xc_func_info_mgga_k_rda;
extern xc_func_info_type xc_func_info_mgga_x_regtm;
extern xc_func_info_type xc_func_info_mgga_k_gea2;
extern xc_func_info_type xc_func_info_mgga_k_gea4;
extern xc_func_info_type xc_func_info_mgga_k_csk1;
extern xc_func_info_type xc_func_info_mgga_k_csk4;
extern xc_func_info_type xc_func_info_mgga_k_csk_loc1;
extern xc_func_info_type xc_func_info_mgga_k_csk_loc4;
extern xc_func_info_type xc_func_info_mgga_k_pc07_opt;
extern xc_func_info_type xc_func_info_mgga_c_kcisk;
extern xc_func_info_type xc_func_info_mgga_c_r2scan01;
extern xc_func_info_type xc_func_info_mgga_c_rmggac;
extern xc_func_info_type xc_func_info_mgga_x_mcml;
extern xc_func_info_type xc_func_info_mgga_x_r2scan01;
extern xc_func_info_type xc_func_info_mgga_x_rppscan;
extern xc_func_info_type xc_func_info_mgga_c_rppscan;
extern xc_func_info_type xc_func_info_mgga_x_r4scan;
extern xc_func_info_type xc_func_info_mgga_x_vcml;
extern xc_func_info_type xc_func_info_mgga_xc_vcml_rvv10;
extern xc_func_info_type xc_func_info_mgga_x_tlda;
extern xc_func_info_type xc_func_info_mgga_x_edmgga;
extern xc_func_info_type xc_func_info_mgga_x_gdme_nv;
extern xc_func_info_type xc_func_info_mgga_x_rlda;
extern xc_func_info_type xc_func_info_mgga_x_gdme_0;
extern xc_func_info_type xc_func_info_mgga_x_gdme_kos;
extern xc_func_info_type xc_func_info_mgga_x_gdme_vt;
extern xc_func_info_type xc_func_info_mgga_x_revtm;
extern xc_func_info_type xc_func_info_mgga_c_revtm;
extern xc_func_info_type xc_func_info_mgga_x_mbrxc_bg;
extern xc_func_info_type xc_func_info_mgga_x_mbrxh_bg;
extern xc_func_info_type xc_func_info_mgga_x_hlta;
extern xc_func_info_type xc_func_info_mgga_c_hltapw;
extern xc_func_info_type xc_func_info_mgga_x_scanl;
extern xc_func_info_type xc_func_info_mgga_x_revscanl;
extern xc_func_info_type xc_func_info_mgga_c_scanl;
extern xc_func_info_type xc_func_info_mgga_c_scanl_rvv10;
extern xc_func_info_type xc_func_info_mgga_c_scanl_vv10;
extern xc_func_info_type xc_func_info_mgga_x_task;
extern xc_func_info_type xc_func_info_mgga_x_mggac;
extern xc_func_info_type xc_func_info_mgga_x_mbr;
extern xc_func_info_type xc_func_info_mgga_x_r2scanl;
extern xc_func_info_type xc_func_info_mgga_c_r2scanl;
extern xc_func_info_type xc_func_info_mgga_x_mtask;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_0;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_1;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_2;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_3;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_4;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_5;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_6;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_7;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_8;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_9;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_10;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_11;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_12;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_13;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_14;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_15;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_16;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_17;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_18;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_19;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_20;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_21;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_22;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_23;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_24;
extern xc_func_info_type xc_func_info_mgga_x_ktbm_gap;
extern xc_func_info_type xc_func_info_mgga_x_mspbel;
extern xc_func_info_type xc_func_info_mgga_x_rmspbel;
extern xc_func_info_type xc_func_info_mgga_x_msrpbel;
extern xc_func_info_type xc_func_info_mgga_x_rmsrpbel;
extern xc_func_info_type xc_func_info_mgga_x_msb86bl;
extern xc_func_info_type xc_func_info_mgga_x_rmsb86bl;

const xc_func_info_type *xc_mgga_known_funct[] = {
  &xc_func_info_mgga_c_dldf,
  &xc_func_info_mgga_xc_zlp,
  &xc_func_info_mgga_xc_otpss_d,
  &xc_func_info_mgga_c_cs,
  &xc_func_info_mgga_c_mn12_sx,
  &xc_func_info_mgga_c_mn12_l,
  &xc_func_info_mgga_c_m11_l,
  &xc_func_info_mgga_c_m11,
  &xc_func_info_mgga_c_m08_so,
  &xc_func_info_mgga_c_m08_hx,
  &xc_func_info_mgga_c_revm11,
  &xc_func_info_mgga_x_lta,
  &xc_func_info_mgga_x_tpss,
  &xc_func_info_mgga_x_m06_l,
  &xc_func_info_mgga_x_gvt4,
  &xc_func_info_mgga_x_tau_hcth,
  &xc_func_info_mgga_x_br89,
  &xc_func_info_mgga_x_bj06,
  &xc_func_info_mgga_x_tb09,
  &xc_func_info_mgga_x_rpp09,
  &xc_func_info_mgga_x_2d_prhg07,
  &xc_func_info_mgga_x_2d_prhg07_prp10,
  &xc_func_info_mgga_x_revtpss,
  &xc_func_info_mgga_x_pkzb,
  &xc_func_info_mgga_x_br89_1,
  &xc_func_info_mgga_k_pgsl025,
  &xc_func_info_mgga_x_ms0,
  &xc_func_info_mgga_x_ms1,
  &xc_func_info_mgga_x_ms2,
  &xc_func_info_mgga_x_th,
  &xc_func_info_mgga_x_m11_l,
  &xc_func_info_mgga_x_mn12_l,
  &xc_func_info_mgga_x_ms2_rev,
  &xc_func_info_mgga_xc_cc06,
  &xc_func_info_mgga_x_mk00,
  &xc_func_info_mgga_c_tpss,
  &xc_func_info_mgga_c_vsxc,
  &xc_func_info_mgga_c_m06_l,
  &xc_func_info_mgga_c_m06_hf,
  &xc_func_info_mgga_c_m06,
  &xc_func_info_mgga_c_m06_2x,
  &xc_func_info_mgga_c_m05,
  &xc_func_info_mgga_c_m05_2x,
  &xc_func_info_mgga_c_pkzb,
  &xc_func_info_mgga_c_bc95,
  &xc_func_info_mgga_c_revtpss,
  &xc_func_info_mgga_xc_tpsslyp1w,
  &xc_func_info_mgga_x_mk00b,
  &xc_func_info_mgga_x_bloc,
  &xc_func_info_mgga_x_modtpss,
  &xc_func_info_mgga_c_tpssloc,
  &xc_func_info_mgga_x_mbeef,
  &xc_func_info_mgga_x_mbeefvdw,
  &xc_func_info_mgga_c_tm,
  &xc_func_info_mgga_xc_b97m_v,
  &xc_func_info_mgga_x_jk,
  &xc_func_info_mgga_x_mvs,
  &xc_func_info_mgga_x_mn15_l,
  &xc_func_info_mgga_c_mn15_l,
  &xc_func_info_mgga_x_scan,
  &xc_func_info_mgga_c_scan,
  &xc_func_info_mgga_c_mn15,
  &xc_func_info_mgga_x_b00,
  &xc_func_info_mgga_xc_hle17,
  &xc_func_info_mgga_c_scan_rvv10,
  &xc_func_info_mgga_x_revm06_l,
  &xc_func_info_mgga_c_revm06_l,
  &xc_func_info_mgga_x_rtpss,
  &xc_func_info_mgga_x_ms2b,
  &xc_func_info_mgga_x_ms2bs,
  &xc_func_info_mgga_x_mvsb,
  &xc_func_info_mgga_x_mvsbs,
  &xc_func_info_mgga_c_revm06,
  &xc_func_info_mgga_c_m06_sx,
  &xc_func_info_mgga_x_ft98,
  &xc_func_info_mgga_c_tpss_gaussian,
  &xc_func_info_mgga_x_eel,
  &xc_func_info_mgga_c_cf22d,
  &xc_func_info_mgga_x_lak,
  &xc_func_info_mgga_c_cc,
  &xc_func_info_mgga_c_ccalda,
  &xc_func_info_mgga_c_rregtm,
  &xc_func_info_mgga_c_b94,
  &xc_func_info_mgga_x_rscan,
  &xc_func_info_mgga_c_rscan,
  &xc_func_info_mgga_x_r2scan,
  &xc_func_info_mgga_c_r2scan,
  &xc_func_info_mgga_x_tm,
  &xc_func_info_mgga_x_vt84,
  &xc_func_info_mgga_x_sa_tpss,
  &xc_func_info_mgga_k_pc07,
  &xc_func_info_mgga_c_kcis,
  &xc_func_info_mgga_xc_lp90,
  &xc_func_info_mgga_c_b88,
  &xc_func_info_mgga_x_gx,
  &xc_func_info_mgga_x_pbe_gx,
  &xc_func_info_mgga_x_revscan,
  &xc_func_info_mgga_c_revscan,
  &xc_func_info_mgga_c_scan_vv10,
  &xc_func_info_mgga_c_revscan_vv10,
  &xc_func_info_mgga_x_br89_explicit,
  &xc_func_info_mgga_x_br89_explicit_1,
  &xc_func_info_mgga_x_regtpss,
  &xc_func_info_mgga_x_2d_js17,
  &xc_func_info_mgga_k_l04,
  &xc_func_info_mgga_k_l06,
  &xc_func_info_mgga_k_rda,
  &xc_func_info_mgga_x_regtm,
  &xc_func_info_mgga_k_gea2,
  &xc_func_info_mgga_k_gea4,
  &xc_func_info_mgga_k_csk1,
  &xc_func_info_mgga_k_csk4,
  &xc_func_info_mgga_k_csk_loc1,
  &xc_func_info_mgga_k_csk_loc4,
  &xc_func_info_mgga_k_pc07_opt,
  &xc_func_info_mgga_c_kcisk,
  &xc_func_info_mgga_c_r2scan01,
  &xc_func_info_mgga_c_rmggac,
  &xc_func_info_mgga_x_mcml,
  &xc_func_info_mgga_x_r2scan01,
  &xc_func_info_mgga_x_rppscan,
  &xc_func_info_mgga_c_rppscan,
  &xc_func_info_mgga_x_r4scan,
  &xc_func_info_mgga_x_vcml,
  &xc_func_info_mgga_xc_vcml_rvv10,
  &xc_func_info_mgga_x_tlda,
  &xc_func_info_mgga_x_edmgga,
  &xc_func_info_mgga_x_gdme_nv,
  &xc_func_info_mgga_x_rlda,
  &xc_func_info_mgga_x_gdme_0,
  &xc_func_info_mgga_x_gdme_kos,
  &xc_func_info_mgga_x_gdme_vt,
  &xc_func_info_mgga_x_revtm,
  &xc_func_info_mgga_c_revtm,
  &xc_func_info_mgga_x_mbrxc_bg,
  &xc_func_info_mgga_x_mbrxh_bg,
  &xc_func_info_mgga_x_hlta,
  &xc_func_info_mgga_c_hltapw,
  &xc_func_info_mgga_x_scanl,
  &xc_func_info_mgga_x_revscanl,
  &xc_func_info_mgga_c_scanl,
  &xc_func_info_mgga_c_scanl_rvv10,
  &xc_func_info_mgga_c_scanl_vv10,
  &xc_func_info_mgga_x_task,
  &xc_func_info_mgga_x_mggac,
  &xc_func_info_mgga_x_mbr,
  &xc_func_info_mgga_x_r2scanl,
  &xc_func_info_mgga_c_r2scanl,
  &xc_func_info_mgga_x_mtask,
  &xc_func_info_mgga_x_ktbm_0,
  &xc_func_info_mgga_x_ktbm_1,
  &xc_func_info_mgga_x_ktbm_2,
  &xc_func_info_mgga_x_ktbm_3,
  &xc_func_info_mgga_x_ktbm_4,
  &xc_func_info_mgga_x_ktbm_5,
  &xc_func_info_mgga_x_ktbm_6,
  &xc_func_info_mgga_x_ktbm_7,
  &xc_func_info_mgga_x_ktbm_8,
  &xc_func_info_mgga_x_ktbm_9,
  &xc_func_info_mgga_x_ktbm_10,
  &xc_func_info_mgga_x_ktbm_11,
  &xc_func_info_mgga_x_ktbm_12,
  &xc_func_info_mgga_x_ktbm_13,
  &xc_func_info_mgga_x_ktbm_14,
  &xc_func_info_mgga_x_ktbm_15,
  &xc_func_info_mgga_x_ktbm_16,
  &xc_func_info_mgga_x_ktbm_17,
  &xc_func_info_mgga_x_ktbm_18,
  &xc_func_info_mgga_x_ktbm_19,
  &xc_func_info_mgga_x_ktbm_20,
  &xc_func_info_mgga_x_ktbm_21,
  &xc_func_info_mgga_x_ktbm_22,
  &xc_func_info_mgga_x_ktbm_23,
  &xc_func_info_mgga_x_ktbm_24,
  &xc_func_info_mgga_x_ktbm_gap,
  &xc_func_info_mgga_x_mspbel,
  &xc_func_info_mgga_x_rmspbel,
  &xc_func_info_mgga_x_msrpbel,
  &xc_func_info_mgga_x_rmsrpbel,
  &xc_func_info_mgga_x_msb86bl,
  &xc_func_info_mgga_x_rmsb86bl,
  NULL
};
