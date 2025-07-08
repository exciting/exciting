/*
 Copyright (C) 2015 Narbe Mardirossian and Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_HYB_MGGA_XC_GAS22   658 /* Google Accelerated Science 22 */

typedef struct {
  double c_x[3], c_ss[5], c_os[5];
} hyb_mgga_xc_gas22_params;

#define N_PAR 16
static const char  *names[N_PAR]  = {
  "_cx00",  "_cx01",  "_cx10",
  "_css00", "_css04", "_css10", "_css20", "_css43",
  "_cos00", "_cos10", "_cos20", "_cos21", "_cos60",
  "_alpha", "_beta", "_omega"};
static const char  *desc[N_PAR]   = {
  "u^00 coefficient for exchange",
  "u^01 coefficient for exchange",
  "u^10 coefficient for exchange",
  "u^01 coefficient for same-spin correlation",
  "u^10 coefficient for same-spin correlation",
  "u^20 coefficient for same-spin correlation",
  "u^06 coefficient for same-spin correlation",
  "u^46 coefficient for same-spin correlation",
  "u^00 coefficient for opposite-spin correlation",
  "u^20 coefficient for opposite-spin correlation",
  "u^60 coefficient for opposite-spin correlation",
  "u^62/3 coefficient for opposite-spin correlation",
  "u^22/3 coefficient for opposite-spin correlation",
  "fraction of HF exchange",
  "fraction of short-range exchange",
  "range-separation constant"
};

static const double par_gas22[N_PAR] = {
  0.862139736374172, 0.936993691972698, 0.317533683085033,
  1.0, -4.10753796482853, -5.24218990333846, -1.76643208454076, 7.5380689617542,
  0.805124374375355, 7.98909430970845, -7.54815900595292, 2.00093961824784, -1.76098915061634,
  1.0, -(1.0 - 0.15), 0.3
};

static void
hyb_mgga_xc_gas22_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(hyb_mgga_xc_gas22_params));
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0);

  p->nlc_b = 6.0;
  p->nlc_C = 0.01;
}

#include "maple2c/mgga_exc/hyb_mgga_xc_gas22.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_xc_gas22 = {
  XC_HYB_MGGA_XC_GAS22,
  XC_EXCHANGE_CORRELATION,
  "Google Accelerated Science 22",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Ma2022_eabq0279, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_NEEDS_TAU | XC_FLAGS_VV10 | MAPLE2C_FLAGS,
  1e-13,
  {N_PAR, names, desc, par_gas22, set_ext_params_cpy_cam},
  hyb_mgga_xc_gas22_init, NULL,
  NULL, NULL, &work_mgga,
};
