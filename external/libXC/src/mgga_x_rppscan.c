/*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_RPPSCAN         648 /* r++SCAN exchange */

typedef struct {
  double c2, d, k1, eta;
} mgga_x_rppscan_params;

#define N_PAR 4
static const char *names[N_PAR] = {"_c2", "_d", "_k1", "_eta"};
static const char *desc[N_PAR] = {"c2 parameter", "d parameter",
  "k1 parameter", "eta parameter"};

static const double par_rppscan[N_PAR] = {0.8, 1.24, 0.065, 0.001};

static void
mgga_x_rppscan_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_rppscan_params));
}

#include "maple2c/mgga_exc/mgga_x_rppscan.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_rppscan = {
  XC_MGGA_X_RPPSCAN,
  XC_EXCHANGE,
  "r++SCAN: rSCAN with uniform density limit and coordinate scaling behavior",
  XC_FAMILY_MGGA,
  {&xc_ref_Furness2022_034109, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-11,
  {N_PAR, names, desc, par_rppscan, set_ext_params_cpy},
  mgga_x_rppscan_init, NULL,
  NULL, NULL, &work_mgga
};
