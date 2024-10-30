/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_MSPBEL          761 /* MS-PBEl, a PBE-like meta-GGA exchange */
#define XC_MGGA_X_RMSPBEL         762 /* regularized MS-PBEl */

typedef struct{
  double kappa, c, b, eta;
} mgga_x_mspbel_params;

static void
mgga_x_mspbel_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_mspbel_params));
}

#define N_PAR 4
static const char  *mspbel_names[N_PAR]  = {"_kappa", "_c", "_b", "_eta"};
static const char  *mspbel_desc[N_PAR]   = {
  "kappa parameter",
  "c parameter",
  "exponent b",
  "eta parameter"
};
static const double mspbel_values[N_PAR]  = {0.804, 0.1036199, 1.0, 0.0};
static const double rmspbel_values[N_PAR]  = {0.804, 0.1036199, 1.0, 0.001};

#include "maple2c/mgga_exc/mgga_x_mspbel.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_mspbel = {
  XC_MGGA_X_MSPBEL,
  XC_EXCHANGE,
  "MS-PBEl, a PBE-like meta-GGA exchange",
  XC_FAMILY_MGGA,
  {&xc_ref_Smeets2019_5395, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, mspbel_names, mspbel_desc, mspbel_values, set_ext_params_cpy},
  mgga_x_mspbel_init, NULL,
  NULL, NULL, &work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_rmspbel = {
  XC_MGGA_X_RMSPBEL,
  XC_EXCHANGE,
  "regularized MS-PBEl",
  XC_FAMILY_MGGA,
  {&xc_ref_Cai2024_8611, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, mspbel_names, mspbel_desc, rmspbel_values, set_ext_params_cpy},
  mgga_x_mspbel_init, NULL,
  NULL, NULL, &work_mgga,
};
