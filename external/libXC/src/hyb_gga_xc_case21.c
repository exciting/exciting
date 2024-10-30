/*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_XC_CASE21    390 /* CASE21 */

typedef struct{
  /* set statically */
  int k; /* order of B splines */
  int Nsp; /* number of B splines */
  double knots[14]; /* knot sequence */

  /* adjustable parameters */
  double cx[10]; /* exchange enhancement */
  double cc[10]; /* correlation enhancement */
  double gammax; /* gamma, exchange */
  double gammac; /* gamma, correlation */
  double ax;     /* fraction of exact exchange */
} hyb_gga_xc_case21_params;

#define N_PAR 23
static const char  *names[N_PAR]      = {
  "_cx0", "_cx1", "_cx2", "_cx3", "_cx4", "_cx5", "_cx6", "_cx7", "_cx8", "_cx9",
  "_cc0", "_cc1", "_cc2", "_cc3", "_cc4", "_cc5", "_cc6", "_cc7", "_cc8", "_cc9",
  "_gammax", "_gammac", "_ax"
};
static const char  *desc[N_PAR]       = {
  "cx0 parameter", "cx1 parameter", "cx2 parameter", "cx3 parameter", "cx4 parameter", "cx5 parameter", "cx6 parameter", "cx7 parameter", "cx8 parameter", "cx9 parameter",
  "cc0 parameter", "cc1 parameter", "cc2 parameter", "cc3 parameter", "cc4 parameter", "cc5 parameter", "cc6 parameter", "cc7 parameter", "cc8 parameter", "cc9 parameter",
  "gammax parameter", "gammac parameter", "ax parameter",
};

static const double case21_values[N_PAR]     = {
  // exchange
  0.889402, 0.997849, 1.11912, 1.24555, 1.35175, 1.4474, 1.54252, 1.63761, 1.73269, 1.82777,
  // correlation
  1.14597, 0.998463, 0.860252, 0.730431, 0.597762, 0.457063, 0.30876, 0.155654, 7.45555E-05, -0.1559416,
  MU_PBE/0.8040, 1.0/0.06672455060314922, 0.25
};

GPU_DEVICE_FUNCTION static double xbspline(double u, int ider, const hyb_gga_xc_case21_params * params) {
  assert(ider<=4);

  double result=0.0;
  double temp[5]; /* dimension ider+1 */
  for(int i=0;i<params->Nsp;i++) {
    xc_bspline(i, params->k, u, ider, params->knots, temp);
    result += params->cx[i]*temp[ider];
  }
  //printf("xbspline %i %e = % e\n",ider,u,result);
  return result;
}

GPU_DEVICE_FUNCTION static double cbspline(double u, int ider, const hyb_gga_xc_case21_params * params) {
  assert(ider<=4);

  double result=0.0;
  double temp[5]; /* dimension ider+1 */
  for(int i=0;i<params->Nsp;i++) {
    xc_bspline(i, params->k, u, ider, params->knots, temp);
    result += params->cc[i]*temp[ider];
  }
  //printf("cbspline %i %e = % e\n",ider,u,result);
  return result;
}

static void
case21_set_ext_params(xc_func_type *p, const double *ext_params)
{
  assert(p != NULL);
  hyb_gga_xc_case21_params * params = (hyb_gga_xc_case21_params *) p->params;

  /* Set internal parameters */
  params->k = 3;
  params->Nsp = 10;
  double qmin = -params->k*1.0/(params->Nsp - params->k);
  double qmax = params->Nsp*1.0/(params->Nsp - params->k);
  int nknots = params->Nsp + params->k + 1;
  double dq  = (qmax - qmin)/(nknots-1);
  for(int k=0; k<nknots; k++) {
    params->knots[k] = qmin+k*dq;
  }

  /* External parameters */
  for(int i=0;i<params->Nsp;i++)
    params->cx[i] = ext_params[i];
  for(int i=0;i<params->Nsp;i++)
    params->cc[i] = ext_params[i + params->Nsp];
  params->gammax = ext_params[2*params->Nsp];
  params->gammac = ext_params[2*params->Nsp+1];
  params->ax = ext_params[2*params->Nsp+2];

  /* Exact exchange */
  set_ext_params_exx(p, ext_params);
}

static void
hyb_gga_xc_case21_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(hyb_gga_xc_case21_params));

  xc_hyb_init_hybrid(p, 0.0);
}

#include "maple2c/gga_exc/hyb_gga_xc_case21.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_case21 = {
  XC_HYB_GGA_XC_CASE21,
  XC_EXCHANGE_CORRELATION,
  "CASE21: Constrained And Smoothed semi-Empirical 2021 functional",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Sparrow2022_6896, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, case21_values, case21_set_ext_params},
  hyb_gga_xc_case21_init, NULL,
  NULL, &work_gga, NULL
};
