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

/************************************************************************
 Implements Perdew, Tao, Staroverov & Scuseria 
   meta-Generalized Gradient Approximation.

  Exchange part
************************************************************************/

#define XC_MGGA_X_TPSS          202 /* Perdew, Tao, Staroverov & Scuseria exchange */

/* some parameters */
static FLOAT b=0.40, c=1.59096, e=1.537, kappa=0.804, mu=0.21951;


/* This is Equation (7) from the paper and its derivatives */
static void 
x_tpss_7(int order, FLOAT p, FLOAT z, 
	 FLOAT *qb, FLOAT *dqbdp, FLOAT *dqbdz)
{
  FLOAT a1, a2, h1, h2;
  FLOAT alpha, dalphadp, dalphadz, dqbdalpha;

  /* Eq. (8) */
  a1    = (1.0/z - 1.0);
  h1    = 5.0/3.0;

  alpha = h1*a1*p;

  /* Eq. (7) */
  a2    = sqrt(1.0 + b*alpha*(alpha-1.0));
  h2    = 9.0/20.0;

  *qb   = h2*(alpha - 1.0)/a2 + 2.0*p/3.0;

  if(order < 1) return;              /* And now the derivatives */

  /* Eq. (8) */
  dalphadp = h1*a1;
  dalphadz = -h1*p/(z*z);

  dqbdalpha = h2*(1.0 + 0.5*b*(alpha - 1.0))/(a2*a2*a2);

  *dqbdp = dqbdalpha*dalphadp + 2.0/3.0;
  *dqbdz = dqbdalpha*dalphadz;
}


/* Equation (10) in all it's glory */
static 
void x_tpss_10(int order, FLOAT p, FLOAT z, 
	       FLOAT *x, FLOAT *dxdp, FLOAT *dxdz)
{
  FLOAT x1, dxdp1, dxdz1;
  FLOAT aux1, z2, p2;
  FLOAT qb, dqbdp, dqbdz;
  
  FLOAT a1, a2, a3, h3, a4, a5, a6, d1;

  /* Equation 7 */
  dqbdp = dqbdz = 0.0;
  x_tpss_7(order, p, z, &qb, &dqbdp, &dqbdz);

  z2   = z*z;
  p2   = p*p;
  aux1 = 10.0/81.0;
  
  /* first we handle the numerator */
  x1    = 0.0;

  a1 = 1.0 + z2;                     /* first term  */
  x1 += (aux1 + c*z2/(a1*a1))*p;

  a2 = 146.0/2025.0;                 /* second term */
  x1 += a2*qb*qb;

  a3 = sqrt(0.5*(9.0*z2/25.0 + p2)); /* third term  */
  h3 = -73.0/405;
  x1 += h3*qb*a3;

  a4 = aux1*aux1/kappa;              /* forth term  */
  x1 += a4*p2;

  a5 = 2.0*sqrt(e)*aux1*9.0/25.0;    /* fifth term  */
  x1 += a5*z2;

  a6 = e*mu;                         /* sixth term  */
  x1 += a6*p*p2;

  d1 = 1.0 + sqrt(e)*p;              /* denominator */
  *x  = x1/(d1*d1);

  if(order < 1) return;              /* the derivatives */

  dxdp1 = 0.0;
  dxdz1 = 0.0;

  dxdp1 += aux1 + c*z2/(a1*a1);              /* first term  */
  dxdz1 += c*2.0*z*(1.0 - z2)*p/(a1*a1*a1);
  
  dxdp1 += 2.0*a2*qb*dqbdp;                  /* second term */
  dxdz1 += 2.0*a2*qb*dqbdz;
  
  dxdp1 += h3*(a3*dqbdp + 0.5*qb*p/a3);      /* third term  */
  dxdz1 += h3*(a3*dqbdz + 0.5*qb*(9.0/25.0)*z/a3);
  
  dxdp1 += a4*2.0*p;                         /* forth term  */

  dxdz1 += a5*2.0*z;                         /* fifth term  */
  
  dxdp1 += a6*3.0*p2;                        /* sixth term  */
  
  *dxdp = (dxdp1*d1 - 2.0*sqrt(e)*x1)/(d1*d1*d1);   /* denominator */
  *dxdz = dxdz1/(d1*d1);
}


static void 
func(const XC(mgga_type) *pt, FLOAT x, FLOAT t, FLOAT u, int order,
     FLOAT *f, FLOAT *vrho0, FLOAT *dfdx, FLOAT *dfdt, FLOAT *dfdu,
     FLOAT *d2fdx2, FLOAT *d2fdxt, FLOAT *d2fdt2)
{
  FLOAT ss, pp, a1, zz;
  FLOAT dxdp, dxdz;
  
  ss = X2S*x;
  pp = ss*ss;

  zz = x*x/(8.0*t);

  /* Eq. 10 */
  x_tpss_10(order, pp, zz, &x, &dxdp, &dxdz);

  /* Eq. (5) */
  a1 = kappa/(kappa + x);
  
  *f = 1.0 + kappa*(1.0 - a1);

  if(order < 1) return;

  *dfdx = a1*a1*(dxdp*2.0*ss*X2S + dxdz*x/(4.0*t));
  *dfdt = a1*a1*dxdz*(-zz/t);
}

#include "work_mgga_x.c"

XC(func_info_type) XC(func_info_mgga_x_tpss) = {
  XC_MGGA_X_TPSS,
  XC_EXCHANGE,
  "Tao, Perdew, Staroverov & Scuseria",
  XC_FAMILY_MGGA,
  "J Tao, JP Perdew, VN Staroverov, and G Scuseria, Phys. Rev. Lett. 91, 146401 (2003)\n"
  "JP Perdew, J Tao, VN Staroverov, and G Scuseria, J. Chem. Phys. 120, 6898 (2004)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  NULL, NULL,
  NULL, NULL,        /* this is not an LDA                   */
  work_mgga_x,
};
