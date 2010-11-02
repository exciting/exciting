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

#define XC_GGA_XC_EDF1        165 /* Empirical functionals from Adamson, Gill, and Pople */
#define XC_GGA_X_OPTPBE_VDW   141 /* PBE reparametrization for vdW */

static void
gga_xc_edf1_init(void *p_)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_GGA_X_B88, XC_GGA_C_LYP};
  static FLOAT funcs_coef[4] = {1.030952 - 10.4017 + 8.44793, 10.4017, -8.44793, 1.0};
  XC(gga_type) *p = (XC(gga_type) *)p_;

  gga_init_mix(p, 4, funcs_id, funcs_coef);  

  XC(gga_x_b88_set_params)(p->func_aux[1], 0.0035, 6.0);
  XC(gga_x_b88_set_params)(p->func_aux[2], 0.0042, 6.0);
  XC(gga_c_lyp_set_params)(p->func_aux[3], 0.055, 0.158, 0.25, 0.3505);
}

const XC(func_info_type) XC(func_info_gga_xc_edf1) = {
  XC_GGA_XC_EDF1,
  XC_EXCHANGE_CORRELATION,
  "EDF1",
  XC_FAMILY_GGA,
  "RD Adamson, PMW Gill, and JA Pople, Chem. Phys. Lett. 284 6 (1998)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  gga_xc_edf1_init, 
  NULL, NULL, NULL
};


static void
gga_x_optpbe_vdw_init(void *p_)
{
  static int   funcs_id  [2] = {XC_GGA_X_PBE, XC_GGA_X_RPBE};
  static FLOAT funcs_coef[2] = {1.0 - 0.054732, 0.054732};
  XC(gga_type) *p = (XC(gga_type) *)p_;

  gga_init_mix(p, 2, funcs_id, funcs_coef);  

  XC(gga_x_pbe_set_params) (p->func_aux[0], 1.04804, 0.175519);
  XC(gga_x_rpbe_set_params)(p->func_aux[1], 1.04804, 0.175519);
}

const XC(func_info_type) XC(func_info_gga_x_optpbe_vdw) = {
  XC_GGA_X_OPTPBE_VDW,
  XC_EXCHANGE,
  "Reparametrized PBE for vdW",
  XC_FAMILY_GGA,
  "J Klimes, DR Bowler, and A Michaelides, J. Phys.: Condens. Matter 22, 022201 (2010)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  gga_x_optpbe_vdw_init, 
  NULL, NULL, NULL
};
