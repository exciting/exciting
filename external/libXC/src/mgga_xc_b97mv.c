/*
 Copyright (C) 2015 Narbe Mardirossian and Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_XC_B97M_V        254 /* Mardirossian and Head-Gordon */

typedef struct {
  double c_x[5], c_ss[5], c_os[5];
} mgga_xc_b97_mv_params;

#define N_PAR 15
static const char  *names[N_PAR]  = {
  "_cx00",  "_cx01",  "_cx02", "_cx10", "_cx11",
  "_css00", "_css02", "_css10", "_css32", "_css42",
  "_cos00", "_cos01", "_cos03", "_cos10", "_cos32"};
static const char  *desc[N_PAR]   = {
  "u^00 coefficient for exchange",
  "u^01 coefficient for exchange",
  "u^02 coefficient for exchange",
  "u^10 coefficient for exchange",
  "u^01 coefficient for exchange",
  "u^00 coefficient for same-spin correlation",
  "u^02 coefficient for same-spin correlation",
  "u^10 coefficient for same-spin correlation",
  "u^32 coefficient for same-spin correlation",
  "u^42 coefficient for same-spin correlation",
  "u^00 coefficient for opposite-spin correlation",
  "u^01 coefficient for opposite-spin correlation",
  "u^03 coefficient for opposite-spin correlation",
  "u^10 coefficient for opposite-spin correlation",
  "u^32 coefficient for opposite-spin correlation",
};

static const double par_b97m_v[N_PAR] = {
   1.000,   1.308,   1.901,   0.416,   3.070,
   1.000,  -1.855,  -5.668, -20.497, -20.364,
   1.000,   1.573,  -6.298,   2.535,  -6.427
};

static void
mgga_xc_b97mv_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_xc_b97_mv_params));
  /* Non-local correlation parameters */
  p->nlc_b = 6.0;
  p->nlc_C = 0.01;
}

#include "maple2c/mgga_exc/mgga_xc_b97mv.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_xc_b97m_v = {
  XC_MGGA_XC_B97M_V,
  XC_EXCHANGE_CORRELATION,
  "B97M-V exchange-correlation functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Mardirossian2015_074111, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_VV10 | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_b97m_v, set_ext_params_cpy},
  mgga_xc_b97mv_init, NULL,
  NULL, NULL, &work_mgga,
};
