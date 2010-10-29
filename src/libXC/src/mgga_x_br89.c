/*
 Copyright (C) 2006-2009 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "util.h"

#define XC_MGGA_X_BR89         206 /* Becke-Roussel 89  */
#define XC_MGGA_X_BJ06         207 /* Becke & Johnson correction to Becke-Roussel 89  */
#define XC_MGGA_X_TB09         208 /* Tran & Blaha correction to Becke & Johnson  */
#define XC_MGGA_X_RPP09        209 /* Rasanen, Pittalis, and Proetto correction to Becke & Johnson  */

typedef struct{
  FLOAT c;
} mgga_x_tb09_params;

static FLOAT br89_gamma = 0.8;


static void mgga_x_tb09_init(void *p_)
{
  XC(mgga_type) *p = (XC(mgga_type) *)p_;
  mgga_x_tb09_params *params;

  assert(p->params == NULL);

  switch(p->info->number){
  case XC_MGGA_X_BR89:  p->func = 0; break;
  case XC_MGGA_X_BJ06:  p->func = 1; break;
  case XC_MGGA_X_TB09:  p->func = 2; break;
  case XC_MGGA_X_RPP09: p->func = 3; break;
  }

  p->params = malloc(sizeof(mgga_x_tb09_params));
  params = (mgga_x_tb09_params *) (p->params);

  /* value of c in Becke-Johnson */
  XC(mgga_x_tb09_set_params_)(p, 1.0);
}


static void mgga_x_tb09_end(void *p_)
{
  XC(mgga_type) *p = (XC(mgga_type) *)p_;

  assert(p->params != NULL);
  free(p->params);
  p->params = NULL;

}


void XC(mgga_x_tb09_set_params)(XC(func_type) *p, FLOAT c)
{
  assert(p != NULL && p->mgga != NULL);
  XC(mgga_x_tb09_set_params_)(p->mgga, c);
}

void XC(mgga_x_tb09_set_params_)(XC(mgga_type) *p, FLOAT c)
{
  mgga_x_tb09_params *params;

  assert(p->params != NULL);
  params = (mgga_x_tb09_params *) (p->params);

  params->c = c;
}

/* This code follows the inversion done in the PINY_MD package */
FLOAT inline 
br_newt_raph(FLOAT a, FLOAT tol,  FLOAT * res, int *ierr)
{
  int count;
  double x, f;
  static int max_iter = 50;

   *ierr = 1;
   if(a == 0.0)
     return 0.0;
   
   /* starting point */
   x = (a < 0.0) ? -1.0 : 1.0;

   count = 0;
   do {
     double arg, eee, xm2, fp;

     xm2 = x - 2.0;
     arg = 2.0*x/3.0;
     eee = exp(-arg)/a;

     f  = x*eee - xm2;
     fp = eee*(1.0 - 2.0/3.0*x) - 1.0;

     x -= f/fp;
     x  = fabs(x);

     count ++;
     *res = fabs(f);
   } while((*res > tol) && (count < max_iter));

   if(count == max_iter) *ierr=0; 
   return x;
}

FLOAT XC(mgga_x_br89_get_x)(FLOAT Q)
{
  FLOAT rhs, br_x, tol, res;
  int ierr;

#if SINGLE_PRECISION
  tol = 1e-6;
#else
  tol = 5e-12;
#endif

  /* build right-hand side of the non-linear equation 
     Remember we use a different definition of tau */
  rhs = 2.0/3.0*POW(M_PI, 2.0/3.0)/Q;

  br_x = br_newt_raph(rhs, tol, &res, &ierr);
  if(ierr == 0){
    fprintf(stderr, 
	    "Warning: Convergence not reached in Becke-Roussel functional\n"
	    "For rhs = %20.14lf (residual = %e)\n", rhs, res);
  }

  return br_x;
}

static void 
func(const XC(mgga_type) *pt, FLOAT x, FLOAT t, FLOAT u, int order,
     FLOAT *f, FLOAT *vrho0, FLOAT *dfdx, FLOAT *dfdt, FLOAT *dfdu,
     FLOAT *d2fdx2, FLOAT *d2fdxt, FLOAT *d2fdt2)
{
  FLOAT Q, br_x, dfdbx, dxdu, ff, dff, v_BR;
  FLOAT cnst, exp1, exp2;

  Q  = (u - 2.0*br89_gamma*t + 0.5*br89_gamma*x*x)/6.0;
  br_x = XC(mgga_x_br89_get_x)(Q);

  cnst = -2.0*POW(M_PI, 1.0/3.0)/X_FACTOR_C;
  exp1 = exp(br_x/3.0);
  exp2 = exp(-br_x);

  v_BR = (ABS(br_x) > MIN_TAU) ?
    exp1*(1.0 - exp2*(1.0 + br_x/2.0))/br_x :
    1.0/2.0 + br_x/6.0 - br_x*br_x/18.0;

  v_BR *= cnst;

  if(pt->func == 0){ /* XC_MGGA_X_BR89 */
    /* we have also to include the factor 1/2 from Eq. (9) */
    *f = - v_BR / 2.0;
  }else{ /* XC_MGGA_X_BJ06 & XC_MGGA_X_TB09 */
    *f = 0.0;
  }

  if(order < 1) return;

  if(pt->func == 0){ /* XC_MGGA_X_BR89 */
    dfdbx = (ABS(br_x) > MIN_TAU) ?
      (3.0 + br_x*(br_x + 2.0) + (br_x - 3.0)/exp2) / (3.0*exp1*exp1*br_x*br_x) :
      1.0/6.0 - br_x/9.0;
    dfdbx *= cnst;
  
    ff  = br_x*exp(-2.0/3.0*br_x)/(br_x - 2);
    dff = -2.0/3.0 + 1.0/br_x - 1.0/(br_x - 2.0); /* dff / ff */

    dxdu  = -1.0/(6.0*Q*dff);

    *dfdx =    x*br89_gamma*dfdbx*dxdu / 2.0;
    *dfdt = -2.0*br89_gamma*dfdbx*dxdu / 2.0;
    *dfdu =                 dfdbx*dxdu / 2.0;

  }else{
    if(pt->func == 1 || pt->func == 2) { /* XC_MGGA_X_BJ0 & XC_MGGA_X_TB09 */
      FLOAT c;
    
      assert(pt->params != NULL);
      c = ((mgga_x_tb09_params *) (pt->params))->c;
    
      *vrho0 = - c*v_BR - (3.0*c - 2.0)*sqrt(5.0/12.0)*sqrt(t)/(X_FACTOR_C*M_PI);

    }else{ /* XC_MGGA_X_RPP09 */
      *vrho0 = - v_BR - sqrt(5.0/12.0)*sqrt(max(t - x*x/4.0, 0.0))/(X_FACTOR_C*M_PI);
    }
  }
}

#include "work_mgga_x.c"

const XC(func_info_type) XC(func_info_mgga_x_br89) = {
  XC_MGGA_X_BR89,
  XC_EXCHANGE,
  "Becke-Roussel 89",
  XC_FAMILY_MGGA,
  "AD Becke and MR Roussel, Phys. Rev. A 39, 3761 (1989)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  NULL, NULL,
  NULL, NULL,        /* this is not an LDA                   */
  work_mgga_x,
};

const XC(func_info_type) XC(func_info_mgga_x_bj06) = {
  XC_MGGA_X_BJ06,
  XC_EXCHANGE,
  "Becke & Johnson 06",
  XC_FAMILY_MGGA,
  "AD Becke and ER Johnson, J. Chem. Phys. 124, 221101 (2006)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_VXC,
  mgga_x_tb09_init,
  mgga_x_tb09_end,
  NULL, NULL,        /* this is not an LDA                   */
  work_mgga_x,
};

const XC(func_info_type) XC(func_info_mgga_x_tb09) = {
  XC_MGGA_X_TB09,
  XC_EXCHANGE,
  "Tran & Blaha 89",
  XC_FAMILY_MGGA,
  "F Tran and P Blaha, Phys. Rev. Lett. 102, 226401 (2009)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_VXC,
  mgga_x_tb09_init,
  mgga_x_tb09_end,
  NULL, NULL,        /* this is not an LDA                   */
  work_mgga_x,
};

const XC(func_info_type) XC(func_info_mgga_x_rpp09) = {
  XC_MGGA_X_RPP09,
  XC_EXCHANGE,
  "Rasanen, Pittalis & Proetto 09",
  XC_FAMILY_MGGA,
  "E Rasanen, S Pittalis & C Proetto, arXiv:0909.1477 (2009)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_VXC,
  mgga_x_tb09_init,
  mgga_x_tb09_end,
  NULL, NULL,        /* this is not an LDA                   */
  work_mgga_x,
};
