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

#define XC_GGA_C_WI  148 /* Wilson & Ivanov */

static inline void 
func(const XC(gga_type) *p, int order, FLOAT rs, FLOAT zeta, FLOAT xt, FLOAT *xs,
     FLOAT *f, FLOAT *dfdrs, FLOAT *dfdz, FLOAT *dfdxt, FLOAT *dfdxs,
     FLOAT *d2fdrs2, FLOAT *d2fdrsz, FLOAT *d2fdrsxt, FLOAT *d2fdrsxs, FLOAT *d2fdz2, 
     FLOAT *d2fdzxt, FLOAT *d2fdzxs, FLOAT *d2fdxt2, FLOAT *d2fdxtxs, FLOAT *d2fdxs2)
{
  const FLOAT a=-0.00652, b=0.0007, c=0.21, d=0.002, k=0.001;

  FLOAT xt0, xt2, xt52, xt72, cnst_rs, num, den;

  cnst_rs = CBRT(4.0*M_PI/3.0);
  xt2     = xt*xt;
  xt0     = sqrt(xt);
  xt52    = xt2*xt0;
  xt72    = xt*xt52;

  num = a + b*xt2*exp(-k*xt2);
  den = c + rs*(1.0 + d*cnst_rs*xt72);

  *f  = num/den;

  if(order < 1) return;

  *dfdrs   = -(1.0 + d*cnst_rs*xt72)*(*f)/den;
  *dfdz    = 0.0;
  *dfdxt   = (-7.0/2.0*d*cnst_rs*rs*xt52*(*f) + 2.0*b*xt*(1.0 - k*xt2)*exp(-k*xt2))/den;
  dfdxs[0] = 0.0;
  dfdxs[1] = 0.0;

  if(order < 2) return;

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
}

#include "work_gga_c.c"

const XC(func_info_type) XC(func_info_gga_c_wi) = {
  XC_GGA_C_WI,
  XC_CORRELATION,
  "Wilson & Ivanov",
  XC_FAMILY_GGA,
  "LC Wilson & S Ivanov, Int. J. Quantum Chem. 69, 523-532 (1998)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  NULL, NULL, NULL,
  work_gga_c
};
