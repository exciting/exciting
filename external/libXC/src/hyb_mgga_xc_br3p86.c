/*
 Copyright (C) 2022 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"
#include "xc_funcs.h"

#define XC_HYB_MGGA_XC_BR3P86       389 /* BR3P86 hybrid meta-GGA from Neumann and Handy */

#define N_PAR 3
static const char  *names[N_PAR]  = {"_a", "_b", "_c"};
static const char  *desc[N_PAR]   = {
  "Fraction of exact exchange",
  "Fraction of BR exchange",
  "Weight for P86 correlation"
};
static const double br3p86_values[N_PAR] = {0.22, 0.67, 0.85};

static void
hyb_mgga_xc_br3p86_init(xc_func_type *p)
{
  /*
    The functional form is not fully described in the Neumann-Handy
    paper; however, they state they use VWN correlation; also the P86
    correction is on top of VWN. According to Aron Cohen, the version
    of VWN in CADPAC is version 5, which is LDA_C_VWN in libxc.
   */

  static int   funcs_id   [4] = {XC_LDA_X, XC_MGGA_X_BR89_1, XC_LDA_C_VWN, XC_GGA_C_P86VWN};
  static double funcs_coef[4] = {0.0, 0.0, 0.0, 0.0}; /* initialized by ext_params */

  xc_mix_init(p, 4, funcs_id, funcs_coef);
  xc_hyb_init_hybrid(p, 0.0);  /* initialized by ext_params */
}

static void
br3p86_set_ext_params(xc_func_type *p, const double *ext_params)
{
  double a, b, c;

  assert(p != NULL);

  a = get_ext_param(p, ext_params, 0);
  b = get_ext_param(p, ext_params, 1);
  c = get_ext_param(p, ext_params, 2);

  p->cam_alpha = a;

  p->mix_coef[0] = 1.0 - a - b;
  p->mix_coef[1] = b;
  p->mix_coef[2] = 1.0 - c;
  p->mix_coef[3] = c;
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_xc_br3p86 = {
  XC_HYB_MGGA_XC_BR3P86,
  XC_EXCHANGE_CORRELATION,
  "BR3P86 hybrid meta-GGA from Neumann and Handy",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Neumann1995_381, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, br3p86_values, br3p86_set_ext_params},
  hyb_mgga_xc_br3p86_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};
