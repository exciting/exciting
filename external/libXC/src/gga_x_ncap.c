/*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"
#include "xc_funcs.h"

#define XC_GGA_X_NCAP  180 /* Nearly correct asymptotic potential */
#define XC_GGA_XC_NCAP 181 /* Nearly correct asymptotic potential + P86 correlation */
#define XC_GGA_X_NCAPR 324 /* Nearly correct asymptotic potential revised */

typedef struct{
  double alpha, beta, mu, zeta;
} gga_x_ncap_params;

static void
gga_x_ncap_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_ncap_params));
}

#include "maple2c/gga_exc/gga_x_ncap.c"
#include "work_gga.c"

#define N_PAR 4
static const char  *names[N_PAR]  = {"_alpha", "_beta", "_mu", "_zeta"};
static const char  *desc[N_PAR]   = {"alpha", "beta", "mu", "zeta"};
/* Precise values of beta, mu and zeta are taken directly from the
   reference implementation in NWChem, while alpha:=4*Pi/3*beta/mu
   according to the paper; it's been evaluated in Maple with 16 digit
   accuracy.
*/
static const double ncap_values[N_PAR] = {
  0.3451117169263783L, 0.01808569669L, 0.2195149727645171L, 0.30412141859531383L
};
/* Precise values of beta, mu and zeta are taken directly from the
   reference implementation obtained from Javier Carmona, while
   alpha:=4*Pi/3*beta/mu according to the paper; it's been evaluated
   in Maple with 16 digit accuracy.
*/
static const double ncapr_values[N_PAR] = {
 0.3431510377915187L, 0.017982946634292535L, 0.2195149727645171L, 0.5L
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ncap = {
  XC_GGA_X_NCAP,
  XC_EXCHANGE,
  "Nearly correct asymptotic potential",
  XC_FAMILY_GGA,
  {&xc_ref_Carmona2019_303, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, ncap_values, set_ext_params_cpy},
  gga_x_ncap_init, NULL,
  NULL, &work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ncapr = {
  XC_GGA_X_NCAPR,
  XC_EXCHANGE,
  "Nearly correct asymptotic potential revised",
  XC_FAMILY_GGA,
  {&xc_ref_Carmona2022_114109, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, ncapr_values, set_ext_params_cpy},
  gga_x_ncap_init, NULL,
  NULL, &work_gga, NULL
};

/* This is how the functional was actually used in the paper */
void
xc_gga_xc_ncap_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_NCAP, XC_GGA_C_P86};
  static double funcs_coef[2] = {1.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_ncap = {
  XC_GGA_XC_NCAP,
  XC_EXCHANGE_CORRELATION,
  "NCAP exchange + P86 correlation",
  XC_FAMILY_GGA,
  {&xc_ref_Carmona2019_303, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  xc_gga_xc_ncap_init, NULL,
  NULL, NULL, NULL
};
