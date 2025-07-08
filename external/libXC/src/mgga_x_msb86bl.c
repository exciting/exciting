/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_MSB86BL          765 /* MS-B86bl, a B86b-like meta-GGA exchange */
#define XC_MGGA_X_RMSB86BL         766 /* regularized MS-B86bl */

typedef struct{
  double kappa, c, b, eta;
} mgga_x_msb86bl_params;

static void
mgga_x_msb86bl_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_msb86bl_params));
}

#define N_PAR 4
static const char  *msb86bl_names[N_PAR]  = {"_kappa", "_c", "_b", "_eta"};
static const char  *msb86bl_desc[N_PAR]   = {
  "kappa parameter",
  "c parameter",
  "exponent b",
  "eta parameter"
};
static const double msb86bl_values[N_PAR]  = {0.804, 0.08809161, 1.0, 0.0};
static const double rmsb86bl_values[N_PAR]  = {0.804, 0.08809161, 1.0, 0.001};

#include "maple2c/mgga_exc/mgga_x_msb86bl.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_msb86bl = {
  XC_MGGA_X_MSB86BL,
  XC_EXCHANGE,
  "MS-B86bl, a B86b-like meta-GGA exchange",
  XC_FAMILY_MGGA,
  {&xc_ref_Smeets2019_5395, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, msb86bl_names, msb86bl_desc, msb86bl_values, set_ext_params_cpy},
  mgga_x_msb86bl_init, NULL,
  NULL, NULL, &work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_rmsb86bl = {
  XC_MGGA_X_RMSB86BL,
  XC_EXCHANGE,
  "regularized MS-B86bl",
  XC_FAMILY_MGGA,
  {&xc_ref_Cai2024_8611, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, msb86bl_names, msb86bl_desc, rmsb86bl_values, set_ext_params_cpy},
  mgga_x_msb86bl_init, NULL,
  NULL, NULL, &work_mgga,
};
