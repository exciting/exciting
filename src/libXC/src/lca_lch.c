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

#define XC_LCA_LCH       302   /* Lee, Colwell & Handy         */

static XC(func_info_type) func_info_lca_lch = {
  XC_LCA_LCH,
  XC_EXCHANGE_CORRELATION,
  "Lee, Colwell & Handy parametrization",
  XC_FAMILY_LCA,
  "A.M. Lee, S.M. Colwell and N.C. Handy, Chem. Phys. Lett. 229, 225 (1994)"
  "A.M. Lee, S.M. Colwell and N.C. Handy, J. Chem. Phys. 103, 10095 (1995)",
  XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC
};


void XC(lca_lch_init)(XC(lca_type) *p)
{
  p->info = &func_info_lca_lch;
}


/* This routine computes the ratio of the orbital susceptibilities fo the interactiog and 
   non-interacting electron gas and its derivative */
void XC(lca_s_lch)(FLOAT rs, FLOAT *s, FLOAT *dsdrs)
{
  static FLOAT c[3] = {1.0, 0.028, -0.042};
  FLOAT tmp;

  tmp    = exp(c[2]*rs);
  *s     = (c[0] + c[1]*rs)*tmp;
  *dsdrs = (c[1] + c[2]*(c[0] + c[1]*rs))*tmp;
}
