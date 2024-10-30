/*
 Copyright (C) 2008 Georg Madsen

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_AIRY  192 /* Constantin et al based on the Airy gas */
#define XC_GGA_X_LAG   193 /* Local Airy Gas */

typedef struct{
  double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
} gga_x_airy_params;

static void
gga_x_airy_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_airy_params));
}

#define N_PAR 10
static const char  *names[N_PAR]  = {"_a1", "_a2", "_a3", "_a4", "_a5", "_a6", "_a7", "_a8", "_a9", "_a10"};
static const char  *desc[N_PAR]   = {
  "a1 parameter",
  "a2 parameter",
  "a3 parameter",
  "a4 parameter",
  "a5 parameter",
  "a6 parameter",
  "a7 parameter",
  "a8 parameter",
  "a9 parameter",
  "a10 parameter"
};

static const double airy_values[N_PAR] =
  {0.041106, 2.626712, 0.092070, 0.657946, 133.983631, 3.217063, 136.707378, 3.223476, 2.675484, 3.473804};
static const double lag_values[N_PAR] =
  {0.041106, 2.626712, 0.092070, 0.657946, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

#include "maple2c/gga_exc/gga_x_airy.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_airy = {
  XC_GGA_X_AIRY,
  XC_EXCHANGE,
  "Constantin et al based on the Airy gas",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2009_035125, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, airy_values, set_ext_params_cpy},
  gga_x_airy_init, NULL,
  NULL, &work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_lag = {
  XC_GGA_X_LAG,
  XC_EXCHANGE,
  "Local Airy Gas",
  XC_FAMILY_GGA,
  {&xc_ref_Vitos2000_10046, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, lag_values, set_ext_params_cpy},
  gga_x_airy_init, NULL,
  NULL, &work_gga, NULL
};
