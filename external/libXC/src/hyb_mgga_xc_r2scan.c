/*
 Copyright (C) 2022 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"
#include "xc_funcs.h"

#define XC_HYB_MGGA_XC_R2SCANH       659 /* r2SCAN hybrid like TPSSh with 10% exact exchange */
#define XC_HYB_MGGA_XC_R2SCAN0       660 /* r2SCAN hybrid like PBE0 with 25% exact exchange */
#define XC_HYB_MGGA_XC_R2SCAN50      661 /* r2SCAN hybrid like PBE50 with 50% exact exchange */


#define N_PAR 1
static const char  *names[N_PAR]      = {"_cx"};
static const char  *desc[N_PAR]       = {"Fraction of exact exchange"};
static const double r2scanh_values[N_PAR]     = {0.10};
static const double r2scan0_values[N_PAR]     = {0.25};
static const double r2scan50_values[N_PAR]    = {0.50};

static void
hyb_mgga_xc_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_MGGA_X_R2SCAN, XC_MGGA_C_R2SCAN};
  static double funcs_coef[2] = {0.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_hybrid(p, 0.0);
}

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  double cx;

  assert(p != NULL);
  cx = get_ext_param(p, ext_params, 0);

  p->mix_coef[0] = 1.0 - cx;
  p->cam_alpha = cx;
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_xc_r2scanh = {
  XC_HYB_MGGA_XC_R2SCANH,
  XC_EXCHANGE_CORRELATION,
  "r2SCANh: r2SCAN hybrid like TPSSh with 10% exact exchange",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Bursch2022_134105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, r2scanh_values, set_ext_params},
  hyb_mgga_xc_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_xc_r2scan0 = {
  XC_HYB_MGGA_XC_R2SCAN0,
  XC_EXCHANGE_CORRELATION,
  "r2SCAN0: r2SCAN hybrid like PBE0 with 25% exact exchange",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Bursch2022_134105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, r2scan0_values, set_ext_params},
  hyb_mgga_xc_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_xc_r2scan50 = {
  XC_HYB_MGGA_XC_R2SCAN50,
  XC_EXCHANGE_CORRELATION,
  "r2SCAN50: r2SCAN hybrid like PBE50 with 50% exact exchange",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Bursch2022_134105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, r2scan50_values, set_ext_params},
  hyb_mgga_xc_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};
