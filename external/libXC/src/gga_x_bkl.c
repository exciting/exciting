/*
 Copyright (C) 2024 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_BKL1        338 /* Bhattacharjee-Koshi-Lee functional for band gaps */
#define XC_GGA_X_BKL2        339 /* Bhattacharjee-Koshi-Lee functional for band gaps */

typedef struct{
  double mu1;
  double kappa;
  double alpha;
  double beta;
  double gamma;
} gga_x_bkl_params;

static void
gga_x_bkl_init(xc_func_type *p)
{
  gga_x_bkl_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_bkl_params));
  params = (gga_x_bkl_params *) (p->params);
}

#define N_PAR 5
static const char  *names[N_PAR]  = {"_mu1", "_kappa", "_alpha", "_beta", "_gamma"};
static const char  *desc[N_PAR]   = {
  "mu1 = mu / kappa",
  "kappa",
  "alpha",
  "beta",
  "gamma"
};

static const double bkl1_values[N_PAR] =
  {MU_PBE/0.8040, 0.8040, 0.01, 2.5, 0.98};
static const double bkl2_values[N_PAR] =
  {MU_PBE/0.8040, 0.8040, 0.01, 0.5, 0.96};

#include "maple2c/gga_exc/gga_x_bkl.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_bkl1 = {
  XC_GGA_X_BKL1,
  XC_EXCHANGE,
  "Type-I band gap functional by Bhattacharjee, Koshi and Lee",
  XC_FAMILY_GGA,
  {&xc_ref_Bhattacharjee2024, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, bkl1_values, set_ext_params_cpy},
  gga_x_bkl_init, NULL,
  NULL, &work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_bkl2 = {
  XC_GGA_X_BKL2,
  XC_EXCHANGE,
  "Type-II band gap functional by Bhattacharjee, Koshi and Lee",
  XC_FAMILY_GGA,
  {&xc_ref_Bhattacharjee2024, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, bkl2_values, set_ext_params_cpy},
  gga_x_bkl_init, NULL,
  NULL, &work_gga, NULL
};
