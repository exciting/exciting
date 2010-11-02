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

#define XC_GGA_X_PBE          101 /* Perdew, Burke & Ernzerhof exchange             */
#define XC_GGA_X_PBE_R        102 /* Perdew, Burke & Ernzerhof exchange (revised)   */
#define XC_GGA_X_PBE_SOL      116 /* Perdew, Burke & Ernzerhof exchange (solids)    */
#define XC_GGA_X_XPBE         123 /* xPBE reparametrization by Xu & Goddard         */
#define XC_GGA_X_PBE_JSJR     126 /* JSJR reparametrization by Pedroza, Silva & Capelle */
#define XC_GGA_X_PBEK1_VDW    140 /* PBE reparametrization for vdW */
#define XC_GGA_X_RGE2         142 /* Regularized PBE */


typedef struct{
  FLOAT kappa, mu;
} gga_x_pbe_params;


static void 
gga_x_pbe_init(void *p_)
{
  static const FLOAT kappa[7] = {
    0.8040,  /* original PBE */
    1.245,   /* PBE R */
    0.8040,  /* PBE sol */
    0.91954, /* xPBE */
    0.8040,  /* PBE_JSJR */
    1.0,     /* PBEK1_VDW */
    0.804    /* RGE2 */
  };

  static const FLOAT mu[7] = {
    0.2195149727645171,   /* PBE: mu = beta*pi^2/3, be ta = 0.066725 */
    0.2195149727645171,   /* PBE rev: as PBE */
    10.0/81.0,            /* PBE sol */
    0.23214,              /* xPBE */
    M_PI*M_PI*0.046/3.0,  /* PBE_JSJR */
    0.2195149727645171,   /* PBEK1_VDW: as PBE */
    10.0/81.0             /* RGE2 */
  };

  XC(gga_type) *p = (XC(gga_type) *)p_;

  assert(p->params == NULL);
  p->params = malloc(sizeof(gga_x_pbe_params));

  switch(p->info->number){
  case XC_GGA_X_PBE_R:      p->func = 1; break;
  case XC_GGA_X_PBE_SOL:    p->func = 2; break;
  case XC_GGA_X_XPBE:       p->func = 3; break;
  case XC_GGA_X_PBE_JSJR:   p->func = 4; break;
  case XC_GGA_X_PBEK1_VDW:  p->func = 5; break;
  case XC_GGA_X_OPTPBE_VDW: p->func = 6; break;
  case XC_GGA_X_RGE2:       p->func = 7; break;
  default:                  p->func = 0; /* original PBE */
  }

  XC(gga_x_pbe_set_params_)(p, kappa[p->func], mu[p->func]);
}


void 
XC(gga_x_pbe_set_params)(XC(func_type) *p, FLOAT kappa, FLOAT mu)
{
  assert(p != NULL && p->gga != NULL);
  XC(gga_x_pbe_set_params_)(p->gga, kappa, mu);
}


void 
XC(gga_x_pbe_set_params_)(XC(gga_type) *p, FLOAT kappa, FLOAT mu)
{
  gga_x_pbe_params *params;

  assert(p->params != NULL);
  params = (gga_x_pbe_params *) (p->params);

  params->kappa = kappa;
  params->mu    = mu;
}


static inline void 
func(const XC(gga_type) *p, int order, FLOAT x, 
     FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  FLOAT kappa, mu, ss, ss2, f0, df0, d2f0;

  assert(p->params != NULL);
  kappa = ((gga_x_pbe_params *) (p->params))->kappa;
  mu    = ((gga_x_pbe_params *) (p->params))->mu;

  ss  = X2S*x;
  ss2 = ss*ss;
 
  f0 = kappa + mu*ss2;
  if(p->info->number == XC_GGA_X_RGE2)
    f0 += mu*mu*ss2*ss2/kappa;

  *f = 1.0 + kappa*(1.0 - kappa/f0);

  if(order < 1) return;

  df0 = 2.0*mu*ss;
  if(p->info->number == XC_GGA_X_RGE2)
    df0 += 4.0*mu*mu*ss2*ss/kappa;

  *dfdx  = X2S*kappa*kappa*df0/(f0*f0);
  *ldfdx = X2S*X2S*mu;

  if(order < 2) return;

  d2f0 = 2.0*mu;
  if(p->info->number == XC_GGA_X_RGE2)
    d2f0 += 4.0*3.0*mu*mu*ss2/kappa;

  *d2fdx2 = X2S*X2S*kappa*kappa/(f0*f0)*(d2f0 - 2.0*df0*df0/f0);
}


void 
XC(gga_x_pbe_enhance)(const XC(gga_type) *p, int order, FLOAT x, 
		      FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  func(p, x, order, f, dfdx, ldfdx, d2fdx2);
}


#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_pbe) = {
  XC_GGA_X_PBE,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof",
  XC_FAMILY_GGA,
  "JP Perdew, K Burke, and M Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)\n"
  "JP Perdew, K Burke, and M Ernzerhof, Phys. Rev. Lett. 78, 1396(E) (1997)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  gga_x_pbe_init, 
  NULL, NULL,
  work_gga_x
};

const XC(func_info_type) XC(func_info_gga_x_pbe_r) = {
  XC_GGA_X_PBE_R,
  XC_EXCHANGE,
  "Revised PBE from Zhang & Yang",
  XC_FAMILY_GGA,
  "Y Zhang and W Yang, Phys. Rev. Lett 80, 890 (1998)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  gga_x_pbe_init, 
  NULL, NULL,
  work_gga_x
};

const XC(func_info_type) XC(func_info_gga_x_pbe_sol) = {
  XC_GGA_X_PBE_SOL,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof SOL",
  XC_FAMILY_GGA,
  "JP Perdew, et al, Phys. Rev. Lett. 100, 136406 (2008)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  gga_x_pbe_init, 
  NULL, NULL,
  work_gga_x
};

const XC(func_info_type) XC(func_info_gga_x_xpbe) = {
  XC_GGA_X_XPBE,
  XC_EXCHANGE,
  "Extended PBE by Xu & Goddard III",
  XC_FAMILY_GGA,
  "X Xu and WA Goddard III, J. Chem. Phys. 121, 4068 (2004)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  gga_x_pbe_init, 
  NULL, NULL,
  work_gga_x
};

const XC(func_info_type) XC(func_info_gga_x_pbe_jsjr) = {
  XC_GGA_X_PBE_JSJR,
  XC_EXCHANGE,
  "Reparametrized PBE by Pedroza, Silva & Capelle",
  XC_FAMILY_GGA,
  "LS Pedroza, AJR da Silva, and K. Capelle, arxiv:0905.1925",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  gga_x_pbe_init, 
  NULL, NULL,
  work_gga_x
};

const XC(func_info_type) XC(func_info_gga_x_pbek1_vdw) = {
  XC_GGA_X_PBEK1_VDW,
  XC_EXCHANGE,
  "Reparametrized PBE for vdW",
  XC_FAMILY_GGA,
  "J Klimes, DR Bowler, and A Michaelides, J. Phys.: Condens. Matter 22, 022201 (2010)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  gga_x_pbe_init, 
  NULL, NULL,
  work_gga_x
};

const XC(func_info_type) XC(func_info_gga_x_rge2) = {
  XC_GGA_X_RGE2,
  XC_EXCHANGE,
  "Regularized PBE",
  XC_FAMILY_GGA,
  "A Ruzsinszky, GI Csonka, and G Scuseria, J. Chem. Theory Comput. 5, 763 (2009)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_x
};

