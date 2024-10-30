/*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_RPPSCAN         649 /* r++SCAN correlation */

typedef struct{
  double eta;   /* regularization parameter */
} mgga_c_rppscan_params;

#define N_PAR 1
static const char  *names[N_PAR]  = {"_eta"};
static const char  *desc[N_PAR]   = {
  "Regularization parameter"};

static const double rppscan_values[N_PAR] = {0.001};
static const double rppscan01_values[N_PAR] = {0.01};
  
#include "maple2c/mgga_exc/mgga_c_rppscan.c"
#include "work_mgga.c"

static void
mgga_c_rppscan_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_c_rppscan_params));
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_rppscan = {
  XC_MGGA_C_RPPSCAN,
  XC_CORRELATION,
  "r++SCAN: rSCAN with uniform density limit and coordinate scaling behavior",
  XC_FAMILY_MGGA,
  {&xc_ref_Furness2022_034109, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, rppscan_values, set_ext_params_cpy},
  mgga_c_rppscan_init, NULL,
  NULL, NULL, &work_mgga,
};
