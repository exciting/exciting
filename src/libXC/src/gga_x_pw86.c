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
#include <assert.h>
#include "util.h"

#define XC_GGA_X_PW86         108 /* Perdew & Wang 86 */

static inline void
func(const XC(gga_type) *p, int order, FLOAT x, 
     FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  const FLOAT aa = 1.296, bb = 14.0, cc = 0.2;
  FLOAT ss, ss2, ss4, dd, d2dd, d3dd;

  ss     = X2S*x;
  ss2    = ss*ss;
  ss4    = ss2*ss2;

  dd     = 1.0 + aa*ss2 + bb*ss4 + cc*ss4*ss2;
  *f     = POW(dd, 1.0/15.0);

  if(order < 1) return;

  d2dd   = ss*(2.0*aa + 4.0*bb*ss2 + 6.0*cc*ss4);

  *dfdx  = X2S*d2dd/15.0 * POW(dd, -14.0/15.0);
  *ldfdx = X2S*X2S*aa/15.0;

  if(order < 2) return;

  d3dd    = 2.0*aa + 4.0*3.0*bb*ss2 + 6.0*5.0*cc*ss4;
  *d2fdx2 = X2S*X2S/15.0 * POW(dd, -14.0/15.0) *
    (d3dd - 14.0/15.0*d2dd*d2dd/dd);
}

#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_pw86) = {
  XC_GGA_X_PW86,
  XC_EXCHANGE,
  "Perdew & Wang 86",
  XC_FAMILY_GGA,
  "JP Perdew and Y Wang, Phys. Rev. B 33, 8800 (1986)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  NULL, NULL, NULL,
  work_gga_x
};
