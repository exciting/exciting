/*
 Copyright (C) 2014 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"
#include "xc_funcs.h"

#define XC_HYB_MGGA_XC_TPSSH       457 /*    TPSS hybrid */
#define XC_HYB_MGGA_XC_REVTPSSH    458 /* revTPSS hybrid */
#define XC_HYB_MGGA_XC_TPSS0       396 /* TPSS hybrid with 25% exact exchange */

#define N_PAR 1
static const char  *names[N_PAR]      = {"_cx"};
static const char  *desc[N_PAR]       = {"Fraction of exact exchange"};
static const double tpssh_values[N_PAR]    = {0.10};
static const double tpss0_values[N_PAR]    = {0.25};

static void
hyb_mgga_xc_tpssh_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_MGGA_X_TPSS, XC_MGGA_C_TPSS};
  static double funcs_coef[2] = {0.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_hybrid(p, 0.0);
}

static void
tpssh_set_ext_params(xc_func_type *p, const double *ext_params)
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
const xc_func_info_type xc_func_info_hyb_mgga_xc_tpssh = {
  XC_HYB_MGGA_XC_TPSSH,
  XC_EXCHANGE_CORRELATION,
  "TPSSh",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Staroverov2003_12129, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, tpssh_values, tpssh_set_ext_params},
  hyb_mgga_xc_tpssh_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};


static void
hyb_mgga_xc_revtpssh_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_MGGA_X_REVTPSS, XC_MGGA_C_REVTPSS};
  static double funcs_coef[2] = {0.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_hybrid(p, 0.0);
}


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_xc_revtpssh = {
  XC_HYB_MGGA_XC_REVTPSSH,
  XC_EXCHANGE_CORRELATION,
  "revTPSSh",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Csonka2010_3688, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, tpssh_values, tpssh_set_ext_params},
  hyb_mgga_xc_revtpssh_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};

static void
hyb_mgga_xc_tpss0_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_MGGA_X_TPSS, XC_MGGA_C_TPSS};
  static double funcs_coef[2] = {0.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_hybrid(p, 0.0);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_xc_tpss0 = {
  XC_HYB_MGGA_XC_TPSS0,
  XC_EXCHANGE_CORRELATION,
  "TPSS0 with 25% exact exchange",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Grimme2005_3067, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, tpss0_values, tpssh_set_ext_params},
  hyb_mgga_xc_tpss0_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};
