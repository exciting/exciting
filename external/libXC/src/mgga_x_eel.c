/*
 Copyright (C) 2024 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_EEL          326 /* Exact Exchange Like response */

typedef struct{
  double c, x0, a0;
} mgga_x_eel_params;

#define N_PAR 3
static const char *names[N_PAR] = {"_c", "_x0", "_a0"};
static const char *desc[N_PAR] = {"c parameter", "x0 parameter", "a0 parameter"};

/* Second parameter is (6*pi)^(-2/3) evaluated with Maple with 20 digits */
static const double par_eel[N_PAR] = {4.759279, 0.14118847627290322299, 3.0};

static void
mgga_x_eel_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_eel_params));
}

#include "maple2c/mgga_exc/mgga_x_eel.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_eel = {
  XC_MGGA_X_EEL,
  XC_EXCHANGE,
  "Exact exchange-like exchange of Aschebrock et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Aschebrock2023_234107, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_eel, set_ext_params_cpy},
  mgga_x_eel_init, NULL,
  NULL, NULL, &work_mgga
};
