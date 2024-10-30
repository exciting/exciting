/*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_FT98     319 /* Filatov and Thiel 1998 */

typedef struct{
  double a;
  double b;
  double a1;
  double a2;
  double b1;
  double b2;
} mgga_x_ft98_params;

#define N_PAR 6

static const char  *names[N_PAR]       = {"_a", "_b", "_a1", "_a2", "_b1", "_b2"};
static const char  *desc[N_PAR]        = {"a", "b", "a1", "a2", "b1", "b2"};
static const double ft98_values[N_PAR] =
  { 0.00528014,  0.00003904539, 2.816049, 0.879058, 0.398773, 66.364138};

static void
mgga_x_ft98_init(xc_func_type *p)
{
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_ft98_params));
}

#include "maple2c/mgga_exc/mgga_x_ft98.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ft98 = {
  XC_MGGA_X_FT98,
  XC_EXCHANGE,
  "Filatov and Thiel 1998 meta-GGA exchange",
  XC_FAMILY_MGGA,
  {&xc_ref_Filatov1998_189, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1.0e-12,
  {N_PAR, names, desc, ft98_values, set_ext_params_cpy},
  mgga_x_ft98_init, NULL,
  NULL, NULL, &work_mgga,
};
