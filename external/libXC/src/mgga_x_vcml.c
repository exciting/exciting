/*
 Copyright (C) 2014 Jess Wellendorff, M.A.L. Marques
               2022 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "xc_funcs.h"

#define XC_MGGA_X_VCML            651 /* VCML exchange, used in VCML-rVV10 by Trepte and Voss */
#define XC_MGGA_XC_VCML_RVV10     652 /* VCML-rVV10 exchange-correlation by Trepte and Voss */

#include "maple2c/mgga_exc/mgga_x_vcml.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_vcml = {
  XC_MGGA_X_VCML,
  XC_EXCHANGE,
  "Exchange part of VCML-rVV10 by Trepte and Voss",
  XC_FAMILY_MGGA,
  {&xc_ref_Trepte2022_1104, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-14,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, &work_mgga,
};

static void
mgga_xc_vcml_rvv10_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_MGGA_X_VCML, XC_GGA_C_REGTPSS};
  static double funcs_coef[2] = {1.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  p->nlc_b = 15.35;
  p->nlc_C = 0.0093;
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_xc_vcml_rvv10 = {
  XC_MGGA_XC_VCML_RVV10,
  XC_EXCHANGE_CORRELATION,
  "VCML-rVV10 by Trepte and Voss",
  XC_FAMILY_MGGA,
  {&xc_ref_Trepte2022_1104, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | XC_FLAGS_I_HAVE_ALL | XC_FLAGS_VV10, /* TBD: this should be rvv10 */
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  mgga_xc_vcml_rvv10_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};
