/*
 Copyright (C) 2022 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_CCALDA           388 /* Iso-orbital corrected LDA correlation */

typedef struct{
  double c; /* parameter in eq 10 */
} mgga_c_ccalda_params;

#define N_PAR 1
static const char  *names[N_PAR]      = {"_c"};
static const char  *desc[N_PAR]       = {"c"};

static const double par_ccalda[N_PAR]    = {10000.0};

#include "maple2c/mgga_exc/mgga_c_ccalda.c"
#include "work_mgga.c"

static void
mgga_c_ccalda_init(xc_func_type *p)
{
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_c_ccalda_params));
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_ccalda = {
  XC_MGGA_C_CCALDA,
  XC_CORRELATION,
  "Iso-orbital corrected LDA correlation by Lebeda et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Lebeda2022_023061, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_ccalda, set_ext_params_cpy},
  mgga_c_ccalda_init, NULL,
  NULL, NULL, &work_mgga
};
