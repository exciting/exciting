/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_X_CAM_S12G         646 /* Swart 2012 range-separated GGA exchange */
#define XC_HYB_GGA_X_CAM_S12H         647 /* Swart 2012 range-separated GGA exchange */

typedef struct {
  double A, B, C, D, E;
} hyb_gga_x_cam_s12_params;

static void
hyb_gga_x_cam_s12_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(hyb_gga_x_cam_s12_params));
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0);
}

#define N_PAR 8
static const char  *names[N_PAR]  = {"_A", "_B", "_C", "_D", "_E", "_alpha", "_beta", "_omega"};
static const char  *desc[N_PAR]   = {
  "A parameter",
  "B parameter",
  "C parameter",
  "D parameter",
  "E parameter",
  "fraction of HF exchange",
  "fraction of SR exchange",
  "range-separation parameter"
};
/* Paper reports alpha and beta in the Yanai convention, which differs from libxc */
static const double par_cam_s12g[N_PAR] = {
  1.03323556, 1.757-1.03323556, 0.00417251, 0.00115216, 0.00706184, 0.34485046, -0.34485046, 1.52420731
};
static const double par_cam_s12h[N_PAR] = {
  1.02149642, 1.757-1.02149642, 0.00825905, 0.00235804, 0.00654977, 0.25+0.10897845, -0.10897845, 0.48516891
};

static void
s12h_set_ext_params(xc_func_type *p, const double *ext_params)
{
  double sr_exx;

  hyb_gga_x_cam_s12_params *params;
  set_ext_params_cpy_cam(p, ext_params);
  params = (hyb_gga_x_cam_s12_params *) (p->params);
}

#include "maple2c/gga_exc/hyb_gga_x_cam_s12.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_x_cam_s12g = {
  XC_HYB_GGA_X_CAM_S12G,
  XC_EXCHANGE,
  "Swart 2012 range-separated hybrid GGA exchange",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Swart2013_166, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_cam_s12g, s12h_set_ext_params},
  hyb_gga_x_cam_s12_init, NULL,
  NULL, &work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_x_cam_s12h = {
  XC_HYB_GGA_X_CAM_S12H,
  XC_EXCHANGE,
  "Swart 2012 range-separated hybrid GGA exchange",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Swart2013_166, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_cam_s12h, s12h_set_ext_params},
  hyb_gga_x_cam_s12_init, NULL,
  NULL, &work_gga, NULL
};
