/*
 Copyright (C) 2006-2007 M.A.L. Marques

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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_B88          106 /* Becke 88 */
#define XC_GGA_X_OPTB88_VDW   139 /* Becke 88 reoptimized to be used with vdW functional of Dion et al*/

typedef struct{
  FLOAT beta, gamma;
} gga_x_b88_params;


static void 
gga_x_b88_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  assert(p->params == NULL);
  p->params = malloc(sizeof(gga_x_b88_params));

  /* value of beta in standard Becke 88 functional */
  switch(p->info->number){
  case XC_GGA_X_OPTB88_VDW:
    p->func = 1; 
    XC(gga_x_b88_set_params_)(p, 0.00336865923905927, 6.98131700797731);
    break;
  default: /* XC_GGA_X_B88 */
    p->func = 0; 
    XC(gga_x_b88_set_params_)(p, 0.0042, 6.0);
    break;
  }
}


static void 
gga_x_b88_end(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  assert(p->params != NULL);
  free(p->params);
  p->params = NULL;
}


void 
XC(gga_x_b88_set_params)(XC(func_type) *p, FLOAT beta, FLOAT gamma)
{
  assert(p != NULL && p->gga != NULL);
  XC(gga_x_b88_set_params_)(p->gga, beta, gamma);
}


void 
XC(gga_x_b88_set_params_)(XC(gga_type) *p, FLOAT beta, FLOAT gamma)
{
  gga_x_b88_params *params;

  assert(p->params != NULL);
  params = (gga_x_b88_params *) (p->params);

  params->beta  = beta;
  params->gamma = gamma;
}


static inline void 
func(const XC(gga_type) *p, int order, FLOAT x, 
     FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  FLOAT f1, f2, df1, df2, d2f1, d2f2;
  FLOAT beta, gamma;

  assert(p->params != NULL);
  beta  = ((gga_x_b88_params *) (p->params))->beta;
  gamma = ((gga_x_b88_params *) (p->params))->gamma;

  f1 = beta/X_FACTOR_C*x*x;
  f2 = 1.0 + gamma*beta*x*asinh(x);
  *f = 1.0 + f1/f2;
 
  if(order < 1) return;

  df1 = 2.0*beta/X_FACTOR_C*x;
  df2 = gamma*beta*(asinh(x) + x/sqrt(1.0 + x*x));

  *dfdx = (df1*f2 - f1*df2)/(f2*f2);
  *ldfdx= beta/X_FACTOR_C;

  if(order < 2) return;

  d2f1 = 2.0*beta/X_FACTOR_C;
  d2f2 = gamma*beta*(2.0 + x*x)/pow(1.0 + x*x, 3.0/2.0);

  *d2fdx2 = (2.0*f1*df2*df2 + d2f1*f2*f2 - f2*(2.0*df1*df2 + f1*d2f2))/(f2*f2*f2);
}

#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_b88) = {
  XC_GGA_X_B88,
  XC_EXCHANGE,
  "Becke 88",
  XC_FAMILY_GGA,
  "AD Becke, Phys. Rev. A 38, 3098 (1988)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  gga_x_b88_init, 
  gga_x_b88_end, 
  NULL,
  work_gga_x
};

const XC(func_info_type) XC(func_info_gga_x_optb88_vdw) = {
  XC_GGA_X_OPTB88_VDW,
  XC_EXCHANGE,
  "opt-Becke 88 for vdW",
  XC_FAMILY_GGA,
  "J Klimes, DR Bowler, and A Michaelides, J. Phys.: Condens. Matter 22, 022201 (2010)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  gga_x_b88_init,
  gga_x_b88_end, 
  NULL,
  work_gga_x
};
