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

#include <assert.h>

#include "util.h"

/* this function converts the spin-density into total density and
	 relative magnetization */
inline void
XC(rho2dzeta)(int nspin, const FLOAT *rho, FLOAT *d, FLOAT *zeta)
{
  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  
  if(nspin==XC_UNPOLARIZED){
    *d    = rho[0];
    *zeta = 0.0;
  }else{
    *d    = rho[0] + rho[1];
    *zeta = (*d > MIN_DENS) ? (rho[0] - rho[1])/(*d) : 0.0;
  }
}
