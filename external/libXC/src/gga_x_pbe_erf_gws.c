/*
 2018 Authored by Andrea Kreppel
 2022 Edited by Henryk Laqua

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 Short-range PBE exchange functional Goll/Werner/Stoll
 E. Goll, H.-J. Werner, and H. Stoll., Phys. Chem. Chem. Phys. 7, 3917 (2005).
 DOI:10.1039/B509242F
*/

#include "util.h"

#define  XC_GGA_X_PBE_ERF_GWS                   655 /* Short ranged PBE exchange (erfc) */
#define  XC_HYB_GGA_X_PBE_ERF_GWS               656 /* Short ranged PBE exchange (erfc) + long-range HF exchange */

#define N_PAR 4

typedef struct{
  double kappa, b_PBE, ax, omega;
} gga_x_pbe_erf_gws_params;


static const char  *param_names[N_PAR]  = {"_kappa", "_b_PBE", "_ax", "_omega"};
static const char  *param_desc[N_PAR]   = {
  "kappa from PBE",
  "original b (AKA mu) from PBE",
  "exponent in short-range b-expansion",
  "range-separation screening parameter (AKA mu)",
};
static const double param_values[N_PAR] = {0.8040, 0.2195149727645171, 19.0, 0.5};


static void xc_gga_x_pbe_erf_gws_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_pbe_erf_gws_params));

  xc_hyb_init_hybrid(p, 0.0);
}

#include "maple2c/gga_exc/gga_x_pbe_erf_gws.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbe_erf_gws = {
  XC_GGA_X_PBE_ERF_GWS,
  XC_EXCHANGE,
  "Short ranged PBE exchange (erfc)",
  XC_FAMILY_GGA,
  {&xc_ref_Goll2005_3917, &xc_ref_Goll2006_276, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, param_names, param_desc, param_values, set_ext_params_cpy_omega},
  xc_gga_x_pbe_erf_gws_init,
  NULL, NULL, &work_gga, NULL
};

static void xc_hyb_gga_x_pbe_erf_gws_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_pbe_erf_gws_params));
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_x_pbe_erf_gws = {
  XC_HYB_GGA_X_PBE_ERF_GWS,
  XC_EXCHANGE,
  "Short-range PBE (GWS) exchange (erfc) + long-range exact exchange",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Goll2005_3917, &xc_ref_Goll2006_276, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, param_names, param_desc, param_values, set_ext_params_cpy_lc},
  xc_hyb_gga_x_pbe_erf_gws_init,
  NULL, NULL, &work_gga, NULL
};
