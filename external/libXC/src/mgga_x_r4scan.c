/*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_R4SCAN         650 /* r4SCAN exchange */

typedef struct {
  double c1, c2, d, k1, eta;
  double dp2, dp4, da4;
} mgga_x_r4scan_params;

#define N_PAR 8
static const char *names[N_PAR] = {"_c1", "_c2", "_d", "_k1", "_eta", "_dp2", "_dp4", "_da4"};
static const char *desc[N_PAR] = {"c1 parameter", "c2 parameter", "d parameter",
  "k1 parameter", "eta parameter", "dp2 parameter", "dp4 parameter", "da4 parameter"};

static const double par_r4scan[N_PAR] = {0.667, 0.8, 1.24, 0.065, 0.001, 0.361, 0.802, 0.178};

static void
mgga_x_r4scan_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_r4scan_params));
}

#include "maple2c/mgga_exc/mgga_x_r4scan.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_r4scan = {
  XC_MGGA_X_R4SCAN,
  XC_EXCHANGE,
  "r$^{4}$SCAN, a functional that satisfies the same exact constraints that SCAN does",
  XC_FAMILY_MGGA,
  {&xc_ref_Furness2022_034109, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-11,
  {N_PAR, names, desc, par_r4scan, set_ext_params_cpy},
  mgga_x_r4scan_init, NULL,
  NULL, NULL, &work_mgga
};
