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

#define XC_GGA_C_AM05         135 /* Armiento & Mattsson 05 correlation             */

static void
gga_c_am05_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_LDA_C_PW_MOD, p->nspin);
}

static void
gga_c_am05_end(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  XC(func_end)(p->func_aux[0]);
  free(p->func_aux[0]);
  free(p->func_aux);
}

static void 
my_gga_c_am05(const void *p_, const FLOAT *rho, const FLOAT *sigma,
	   FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	   FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  const FLOAT am05_alpha = 2.804;
  const FLOAT am05_gamma = 0.8098;

  FLOAT sfact, dens, m_zk, vrho_LDA[2];
  int is;

  XC(gga_type) *p = (XC(gga_type) *)p_;

  sfact = (p->nspin == XC_POLARIZED) ? 1.0 : 2.0;

  dens = rho[0];
  if(p->nspin == XC_POLARIZED) dens += rho[1];
  if(dens < MIN_DENS) return;

  XC(lda_exc_vxc)(p->func_aux[0], 1, rho, &m_zk, vrho_LDA);

  for(is=0; is<p->nspin; is++){
    FLOAT gdm, ds, rho13;
    FLOAT x=0.0, ss=0.0, XX, f;
    int js = is==0 ? 0 : 2;

    gdm   = SQRT(sigma[js])/sfact;  
    ds    = rho[is]/sfact;
    rho13 = CBRT(ds);

    if(rho[is] < MIN_DENS){
      /* This is how I think it should be */
      XX = (vsigma[js] < MIN_GRAD) ? 1.0 : 0.0;
    }else{
      x  = gdm/(ds*rho13);
      ss = x*X2S;
    
      XX = 1.0/(1.0 + am05_alpha*ss*ss);
    }

    f = XX + (1.0 - XX)*am05_gamma;

    if(zk != NULL)
      *zk += sfact*m_zk*ds*f;

    if(vrho != NULL){
      int jj, n;
      FLOAT dXX, df;

      if(rho[is] >= MIN_DENS){
	dXX = -2.0*am05_alpha*ss * XX*XX;
      }else{
	dXX = 0.0;
      }
      df  = dXX*(1.0 - am05_gamma)*X2S;

      n = (p->nspin == XC_POLARIZED) ? 2 : 1;
      for(jj=0; jj<n; jj++)
	vrho[jj] += sfact*(vrho_LDA[jj] - m_zk)/dens * ds*f;
      vrho[is] += m_zk*f;
      if(rho[is] >= MIN_DENS)
	 vrho[is] += -4.0/3.0*m_zk*df*x;

      if(rho[is] >= MIN_DENS && gdm>MIN_GRAD)
	vsigma[js] = sfact*m_zk*ds*df*x/(2.0*sigma[js]);
    }
  }

  if(zk != NULL)
    *zk /= dens;  /* we want energy per particle */
  
}

/* Warning: this is a workaround to support blocks while waiting for the next interface */
static void 
gga_c_am05(const void *p_, int np, const FLOAT *rho, const FLOAT *sigma,
	  FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	  FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  int ip;
  const XC(gga_type) *p = p_;

  for(ip=0; ip<np; ip++){
    my_gga_c_am05(p_, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);

    /* increment pointers */
    rho   += p->n_rho;
    sigma += p->n_sigma;
    
    if(zk != NULL)
      zk += p->n_zk;
    
    if(vrho != NULL){
      vrho   += p->n_vrho;
      vsigma += p->n_vsigma;
    }

    if(v2rho2 != NULL){
      v2rho2     += p->n_v2rho2;
      v2rhosigma += p->n_v2rhosigma;
      v2sigma2   += p->n_v2sigma2;
    }
  }
}

const XC(func_info_type) XC(func_info_gga_c_am05) = {
  XC_GGA_C_AM05,
  XC_CORRELATION,
  "Armiento & Mattsson 05",
  XC_FAMILY_GGA,
  "R Armiento and AE Mattsson, Phys. Rev. B 72, 085108 (2005)\n"
  "AE Mattsson, R Armiento, J Paier, G Kresse, JM Wills, and TR Mattsson, J. Chem. Phys. 128, 084714 (2008).",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  gga_c_am05_init,
  gga_c_am05_end,
  NULL,            /* this is not an LDA                   */
  gga_c_am05,
};
