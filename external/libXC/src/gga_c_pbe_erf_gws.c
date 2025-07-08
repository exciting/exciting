/*
 2018 Authored by Andrea Kreppel
 2022 Edited by Henryk Laqua

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 Short-range PBE correlation functional Goll/Werner/Stoll
 Goll, Werner, Stoll Phys. Chem. Chem. Phys. 7, (2005) 3917.
*/

#include "util.h"

#define  XC_GGA_C_PBE_ERF_GWS                   657 /* Short ranged PBE correlation (erfc) */

#define N_PAR 4

typedef struct{
  double beta, gamma, a_c, omega;
} gga_c_pbe_erf_gws_params;

static void
xc_gga_c_pbe_erf_gws_init(xc_func_type *p)
{
  xc_hyb_init_hybrid(p, 0.0);
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_c_pbe_erf_gws_params));
}

static const char  *param_names[N_PAR]  = {"_beta", "_gamma", "_a_c","_omega"};
static const char  *param_desc[N_PAR]   = {
  "beta constant",
  "(1 - ln(2))/Pi^2 in PBE",
  "exponent in beta expansion",
  "range-separation screening parameter (AKA mu)",
};
static const double param_values[N_PAR] =
  {0.06672455060314922, 0.031090690869654895034, 2.78, 0.5};

#include "maple2c/gga_exc/gga_c_pbe_erf_gws.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_pbe_erf_gws = {
  XC_GGA_C_PBE_ERF_GWS,
  XC_CORRELATION,
  "Short ranged PBE correlation (erfc)",
  XC_FAMILY_GGA,
  {&xc_ref_Goll2005_3917,&xc_ref_Goll2006_276,NULL,NULL,NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-13,
  {N_PAR, param_names, param_desc, param_values, set_ext_params_cpy_omega},
  xc_gga_c_pbe_erf_gws_init, NULL,
  NULL, &work_gga, NULL
};
