/*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_X_YUKAWA   641   /* Short-range LDA exchange with Yukawa attenuation */

#include "maple2c/lda_exc/lda_x_yukawa.c"
#include "work_lda.c"

static void
xc_lda_x_yukawa_init(xc_func_type *p)
{
  xc_hyb_init_hybrid(p, 0.0);
}

static const char  *omega_names[]  = {"_omega"};
static const char  *omega_desc[]   = {"screening parameter"};
static const double omega_values[] = {0.3};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_x_yukawa = {
  XC_LDA_X_YUKAWA,
  XC_EXCHANGE,
  "Short-range LDA exchange with Yukawa attenuation",
  XC_FAMILY_LDA,
  {&xc_ref_Savin1995_327, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-13,
  {1, omega_names, omega_desc, omega_values, set_ext_params_cpy_omega},
  xc_lda_x_yukawa_init, NULL,
  &work_lda, NULL, NULL
};
