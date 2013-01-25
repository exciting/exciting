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

#define XC_GGA_C_LYP  131 /* Lee, Yang & Parr */

typedef struct{
  FLOAT A, B, c, d;
} gga_c_lyp_params;

void XC(gga_c_lyp_init)(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;
  gga_c_lyp_params *params;

  assert(p->params == NULL);

  p->params = malloc(sizeof(gga_c_lyp_params));
  params = (gga_c_lyp_params *) (p->params);

  /* values of constants in standard LYP functional */
  XC(gga_c_lyp_set_params_)(p, 0.04918, 0.132, 0.2533, 0.349);
}


void XC(gga_c_lyp_set_params)(XC(func_type) *p, FLOAT A, FLOAT B, FLOAT c, FLOAT d)
{
  assert(p != NULL && p->gga != NULL);
  XC(gga_c_lyp_set_params_)(p->gga, A, B, c, d);
}


void XC(gga_c_lyp_set_params_)(XC(gga_type) *p, FLOAT A, FLOAT B, FLOAT c, FLOAT d)
{
  gga_c_lyp_params *params;

  assert(p->params != NULL);
  params = (gga_c_lyp_params *) (p->params);

  params->A = A;
  params->B = B;
  params->c = c;
  params->d = d;
}


static inline void 
func(const XC(gga_type) *p, int order, FLOAT rs, FLOAT zeta, FLOAT xt, FLOAT *xs,
     FLOAT *f, FLOAT *dfdrs, FLOAT *dfdz, FLOAT *dfdxt, FLOAT *dfdxs,
     FLOAT *d2fdrs2, FLOAT *d2fdrsz, FLOAT *d2fdrsxt, FLOAT *d2fdrsxs, FLOAT *d2fdz2, 
     FLOAT *d2fdzxt, FLOAT *d2fdzxs, FLOAT *d2fdxt2, FLOAT *d2fdxtxs, FLOAT *d2fdxs2)
{
  gga_c_lyp_params *params;

  FLOAT AA, BB, cc, dd, Cf; /* sortcuts for parameters */
  FLOAT cnst_rs, xt2, opz, omz, opz53, omz53, opz83, omz83, zeta2, opdrs;
  FLOAT omega, delta, domega, ddelta;
  FLOAT aux4, aux5, aux6;
  FLOAT t1, t2, t3, t4, t5, t6;
  FLOAT dt1drs, dt2drs, dt4drs, dt5drs;
  FLOAT dt1dz, dt2dz, dt3dz, dt4dz, dt5dz, dt6dz;

  assert(p->params != NULL);
  params = (gga_c_lyp_params *)(p->params);

  cnst_rs = CBRT(4.0*M_PI/3.0);
  Cf      = 3.0*POW(3.0*M_PI*M_PI, 2.0/3.0)/10.0;
  xt2     = xt*xt;

  AA = params->A;
  BB = params->B;
  cc = params->c*cnst_rs;
  dd = params->d*cnst_rs;

  zeta2 = zeta*zeta;
  opz = 1.0 + zeta;
  omz = 1.0 - zeta;
  opz53 = POW(opz, 5.0/3.0);
  omz53 = POW(omz, 5.0/3.0);   
  opz83 = opz*opz53;
  omz83 = omz*omz53;

  opdrs = 1.0/(1.0 + dd*rs);
  omega = BB*exp(-cc*rs)*opdrs;
  delta = (cc + dd*opdrs)*rs;

  aux6 = 1.0/POW(2.0, 8.0/3.0);
  aux4 = aux6/4.0;
  aux5 = aux4/(9.0*2.0);

  t1 = -(1.0 - zeta2)/(1.0 + dd*rs);
  t2 = -xt2*((1.0 - zeta2)*(47.0 - 7.0*delta)/(4.0*18.0) - 2.0/3.0);
  t3 = -Cf/2.0*(1.0 - zeta2)*(opz83 + omz83);
  t4 =  aux4*(1.0 - zeta2)*(5.0/2.0 - delta/18.0)*(xs[0]*xs[0]*opz83 + xs[1]*xs[1]*omz83);
  t5 =  aux5*(1.0 - zeta2)*(delta - 11.0)*(xs[0]*xs[0]*opz*opz83 + xs[1]*xs[1]*omz*omz83);
  t6 = -aux6*(2.0/3.0*(xs[0]*xs[0]*opz83 + xs[1]*xs[1]*omz83) -
	      opz*opz*xs[1]*xs[1]*omz83/4.0 - omz*omz*xs[0]*xs[0]*opz83/4.0);

  *f  = AA*(t1 + omega*(t2 + t3 + t4 + t5 + t6));

  if(order < 1) return;

  domega = -omega*(cc + dd*opdrs);
  ddelta = cc + dd*opdrs*(1.0 - dd*rs*opdrs);

  dt1drs = -dd*t1/(1.0 + dd*rs);
  dt2drs =  xt2*(1.0 - zeta2)*7.0*ddelta/(4.0*18.0);
  dt4drs = -aux4*(1.0 - zeta2)*ddelta/18.0*(xs[0]*xs[0]*opz83 + xs[1]*xs[1]*omz83);
  dt5drs =  aux5*(1.0 - zeta2)*ddelta*(xs[0]*xs[0]*opz*opz83 + xs[1]*xs[1]*omz*omz83);
  *dfdrs =  AA*(dt1drs + domega*(t2 + t3 + t4 + t5 + t6) + omega*(dt2drs + dt4drs + dt5drs));
  
  dt1dz =  2.0*zeta/(1.0 + dd*rs);
  dt2dz =  xt2*2.0*zeta*(47.0 - 7.0*delta)/(4.0*18.0);
  dt3dz = -Cf/2.0*(-2.0*zeta*(opz83 + omz83) + (1.0 - zeta2)*8.0/3.0*(opz53 - omz53));
  dt4dz =  aux4*(5.0/2.0 - delta/18.0)*
    (-2.0*zeta*(xs[0]*xs[0]*opz83 + xs[1]*xs[1]*omz83) + (1.0 - zeta2)*8.0/3.0*(xs[0]*xs[0]*opz53 - xs[1]*xs[1]*omz53));
  dt5dz =  aux5*(delta - 11.0)*
    (-2.0*zeta*(xs[0]*xs[0]*opz*opz83 + xs[1]*xs[1]*omz*omz83) + (1.0 - zeta2)*11/3.0*(xs[0]*xs[0]*opz83 - xs[1]*xs[1]*omz83));
  dt6dz = -aux6*(2.0/3.0*8.0/3.0*(xs[0]*xs[0]*opz53 - xs[1]*xs[1]*omz53)
		 - xs[1]*xs[1]/4.0*opz*( 2.0*omz83 - opz*8.0/3.0*omz53)
		 - xs[0]*xs[0]/4.0*omz*(-2.0*opz83 + omz*8.0/3.0*opz53));
  *dfdz =  AA*(dt1dz + omega*(dt2dz + dt3dz + dt4dz + dt5dz + dt6dz));

  *dfdxt = -2.0*AA*omega*xt*((1.0 - zeta2)*(47.0 - 7.0*delta)/(4.0*18.0) - 2.0/3.0);

  dfdxs[0] = AA*omega*2.0*xs[0]*
    (aux4*(1.0 - zeta2)*(5.0/2.0 - delta/18.0)*opz83 +
     aux5*(1.0 - zeta2)*(delta - 11.0)*opz*opz83 -
     aux6*(2.0/3.0*opz83 - omz*omz*opz83/4.0)
     );
  dfdxs[1] = AA*omega*2*xs[1]*
    (aux4*(1.0 - zeta2)*(5.0/2.0 - delta/18.0)*omz83 +
     aux5*(1.0 - zeta2)*(delta - 11.0)*omz*omz83 -
     aux6*(2.0/3.0*omz83 - opz*opz*omz83/4.0)
     );

  if(order < 2) return;

  /*
  *d2fdrs2    = -2.0*(1.0 + d*cnst_rs*xt72)*(*dfdrs)/den;
  *d2fdrsz    = 0.0;
  *d2fdrsxt   = -7.0*d*cnst_rs*rs*xt52*(*dfdrs)/den - 7.0/2.0*d*cnst_rs*xt52*(*f)/den
    - 2.0*b*xt*exp(-k*xt2)*(1.0 - k*xt2)*(1.0 + d*cnst_rs*xt72)/(den*den);
  d2fdrsxs[0] = 0.0;
  d2fdrsxs[1] = 0.0;
  *d2fdz2     = 0.0;
  *d2fdzxt    = 0.0;
  d2fdzxs[0]  = 0.0;
  d2fdzxs[1]  = 0.0;
  *d2fdxt2    = (49.0/2.0*d*cnst_rs*rs*xt2*xt2/den - 35.0/4.0*xt0)*(*f)*d*cnst_rs*rs*xt/den
    + (-7.0*d*cnst_rs*rs*xt72*(1.0 - k*xt2)/den + 1.0 - 5.0*k*xt2 + 2.0*k*k*xt2*xt2)*2.0*b*exp(-k*xt2)/den;
  d2fdxtxs[0] = 0.0;
  d2fdxtxs[1] = 0.0;
  d2fdxs2[0]  = 0.0;
  d2fdxs2[1]  = 0.0;
  d2fdxs2[2]  = 0.0;
  */
}

#include "work_gga_c.c"

const XC(func_info_type) XC(func_info_gga_c_lyp) = {
  XC_GGA_C_LYP,
  XC_CORRELATION,
  "Lee, Yang & Parr",
  XC_FAMILY_GGA,
  "C Lee, W Yang and RG Parr, Phys. Rev. B 37, 785 (1988)\n"
  "B Miehlich, A Savin, H Stoll and H Preuss, Chem. Phys. Lett. 157, 200 (1989)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC, // | XC_FLAGS_HAVE_FXC,
  XC(gga_c_lyp_init), 
  NULL,
  NULL,
  work_gga_c
};
