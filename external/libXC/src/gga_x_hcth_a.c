/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_HCTH_A          34 /* HCTH-A */

typedef struct{
  double beta, gamma;
  double c0, c1, c2;
} gga_x_hcth_a_params;


static void
gga_x_hcth_a_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_hcth_a_params));
}

#define N_PAR 5
static const char  *names[N_PAR]  = {"_beta", "_gamma", "_c0", "_c1", "_c2"};
static const char  *desc[N_PAR]   = {
  "beta/X_FACTOR_C is the coefficient of the gradient expansion",
  "gamma should be 6 to get the right asymptotics of Ex",
  "LDA coefficient",
  "B88 coefficient",
  "dB88/dbeta coefficient"
};
static const double hcth_a_values[N_PAR] =
  {0.0042, 6.0, 1.09878, -2.51173, 0.0156233};

#include "maple2c/gga_exc/gga_x_hcth_a.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_hcth_a = {
  XC_GGA_X_HCTH_A,
  XC_EXCHANGE,
  "HCTH-A",
  XC_FAMILY_GGA,
  {&xc_ref_Hamprecht1998_6264, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, hcth_a_values, set_ext_params_cpy},
  gga_x_hcth_a_init, NULL,
  NULL, &work_gga, NULL
};
