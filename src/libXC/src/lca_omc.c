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

#define XC_LCA_OMC       301   /* Orestes, Marcasso & Capelle  */

static XC(func_info_type) func_info_lca_omc = {
  XC_LCA_OMC,
  XC_EXCHANGE_CORRELATION,
  "Orestes, Marcasso & Capelle parametrization",
  XC_FAMILY_LCA,
  "E. Orestes, T. Marcasso, and K. Capelle, Phys. Rev. A 68, 022105 (2003)",
  XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC
};

void XC(lca_omc_init)(XC(lca_type) *p)
{
  p->info = &func_info_lca_omc;
}


/* This routine computes the ratio of the orbital susceptibilities for the interacting and 
   non-interacting electron gas and its derivative */
void XC(lca_s_omc)(FLOAT rs, FLOAT *s, FLOAT *dsdrs)
{
  static FLOAT c[5] = {1.1038, -0.4990, 0.4423, -0.06696, 0.0008432};
  FLOAT tmp;
  
  tmp    = SQRT(rs);
  *s     = c[0] + c[1]*CBRT(rs) + c[2]*tmp + c[3]*rs + c[4]*rs*rs;
  *dsdrs = c[1]*POW(rs, -2.0/3.0)/3.0 + c[2]/(2.0*tmp) + c[3] + 2.0*c[4]*rs;
}
