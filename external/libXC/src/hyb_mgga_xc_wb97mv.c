/*
 Copyright (C) 2015 Narbe Mardirossian and Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_HYB_MGGA_XC_WB97M_V   531 /* Mardirossian and Head-Gordon */

typedef struct {
  double c_x[3], c_ss[5], c_os[6];
} hyb_mgga_xc_wb97_mv_params;

#define N_PAR 17
static const char  *names[N_PAR]  = {
  "_cx00",  "_cx01",  "_cx10",
  "_css00", "_css04", "_css10", "_css20", "_css43",
  "_cos00", "_cos10", "_cos20", "_cos21", "_cos60", "_cos61", 
  "_alpha", "_beta", "_omega"};
static const char  *desc[N_PAR]   = {
  "u^00 coefficient for exchange",
  "u^01 coefficient for exchange",
  "u^10 coefficient for exchange",
  "u^00 coefficient for same-spin correlation",
  "u^04 coefficient for same-spin correlation",
  "u^10 coefficient for same-spin correlation",
  "u^20 coefficient for same-spin correlation",
  "u^43 coefficient for same-spin correlation",
  "u^00 coefficient for opposite-spin correlation",
  "u^10 coefficient for opposite-spin correlation",
  "u^20 coefficient for opposite-spin correlation",
  "u^21 coefficient for opposite-spin correlation",
  "u^60 coefficient for opposite-spin correlation",
  "u^61 coefficient for opposite-spin correlation",
  "fraction of HF exchange",
  "fraction of short-range exchange",
  "range-separation constant"
};

static const double par_wb97m_v[N_PAR] = {
   0.85,   1.007,  0.259,
   0.443, -1.437, -4.535,    -3.39,          4.278,
   1.0,         1.358,        2.924,       -8.812,     -1.39,         9.142, 
   1.0, -(1.0 - 0.15), 0.3
};

static void
hyb_mgga_xc_wb97mv_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(hyb_mgga_xc_wb97_mv_params));
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0);

  p->nlc_b = 6.0;
  p->nlc_C = 0.01;
}

#include "maple2c/mgga_exc/hyb_mgga_xc_wb97mv.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_xc_wb97m_v = {
  XC_HYB_MGGA_XC_WB97M_V,
  XC_EXCHANGE_CORRELATION,
  "wB97M-V exchange-correlation functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Mardirossian2016_214110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_NEEDS_TAU | XC_FLAGS_VV10 | MAPLE2C_FLAGS,
  1e-13,
  {N_PAR, names, desc, par_wb97m_v, set_ext_params_cpy_cam},
  hyb_mgga_xc_wb97mv_init, NULL,
  NULL, NULL, &work_mgga,
};
