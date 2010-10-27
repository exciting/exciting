/*
 Copyright (C) 2006-2008 M.A.L. Marques

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

/************************************************************************
  This file is to be included in meta GGA exchange functionals. As often these
  functionals are written as a function of s = |grad n|/n^(4/3) and tau, this
  routine performs the necessary conversions between a functional of s and tau
  and of rho.
************************************************************************/

#include <stdio.h>

static void 
work_mgga_x(const void *p_, int np,
	    const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
	    FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau,
	    FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2)
{
  const XC(mgga_type) *p = p_;

  FLOAT sfact, sfact2, dens;
  int is, ip, order;

  order = -1;
  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;
  if(order < 0) return;

  sfact = (p->nspin == XC_POLARIZED) ? 1.0 : 2.0;
  sfact2 = sfact*sfact;

  for(ip = 0; ip < np; ip++){
    dens = (p->nspin == XC_UNPOLARIZED) ? rho[0] : rho[0] + rho[1];
    if(dens < MIN_DENS) goto end_ip_loop;

    for(is=0; is<p->nspin; is++){
      FLOAT gdm, ds, rho13;
      FLOAT x, t, u, f, lnr2, ltau, vrho0, dfdx, dfdt, dfdu, d2fdx2, d2fdxt, d2fdt2;
      int js = (is == 0) ? 0 : 2;

      if(rho[is] < MIN_DENS || tau[is] < MIN_TAU) continue;

      gdm   = sqrt(sigma[js])/sfact;
      ds    = rho[is]/sfact;
      rho13 = POW(ds, 1.0/3.0);
      x     = gdm/(ds*rho13);
    
      ltau  = tau[is]/sfact;
      t     = ltau/(ds*rho13*rho13);  /* tau/rho^(5/3) */

      lnr2  = lapl_rho[is]/sfact;     /* this can be negative */
      u     = lnr2/(ds*rho13*rho13);  /* lapl_rho/rho^(5/3) */

      vrho0 = dfdx = dfdt = dfdu = 0.0;
      d2fdx2 = d2fdxt = d2fdt2 = 0.0;

      func(p, x, t, u, order, &f, &vrho0,
	   &dfdx, &dfdt, &dfdu, &d2fdx2, &d2fdxt, &d2fdt2);

      if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
	*zk += -sfact*X_FACTOR_C*(ds*rho13)*f;

      if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC)){
	vrho[is]      = -X_FACTOR_C*rho13*(vrho0 + 4.0/3.0*(f - dfdx*x) - 5.0/3.0*(dfdt*t + dfdu*u));
	vtau[is]      = -X_FACTOR_C*dfdt/rho13;
	vlapl_rho[is] = -X_FACTOR_C*dfdu/rho13;
	if(gdm>MIN_GRAD)
	  vsigma[js]    = -sfact*X_FACTOR_C*(rho13*ds)*dfdx*x/(2.0*sigma[js]);
      }

      if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC)){
	/* Missing terms here */
	exit(1);
      }
    }
    
    if(zk != NULL)
      *zk /= dens; /* we want energy per particle */

  end_ip_loop:
    /* increment pointers */
    rho      += p->n_rho;
    sigma    += p->n_sigma;
    tau      += p->n_tau;
    lapl_rho += p->n_lapl_rho;
    
    if(zk != NULL)
      zk += p->n_zk;
    
    if(vrho != NULL){
      vrho      += p->n_vrho;
      vsigma    += p->n_vsigma;
      vtau      += p->n_vtau;
      vlapl_rho += p->n_vlapl_rho;
    }

    if(v2rho2 != NULL){
      v2rho2     += p->n_v2rho2;
      v2rhosigma += p->n_v2rhosigma;
      v2sigma2   += p->n_v2sigma2;
      /* warning: extra termns missing */
    }
  }
}
