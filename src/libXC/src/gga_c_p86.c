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

#define XC_GGA_C_P86 132 /* Perdew 86 */

/************************************************************************
 Implements Perdew 86 Generalized Gradient Approximation
 correlation functional.
************************************************************************/

/* TODO: convert to perdew functionals */

static void
gga_c_p86_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  p->n_func_aux  = 1;
  p->func_aux    = (XC(func_type) **) malloc(sizeof(XC(func_type) *)*p->n_func_aux);
  p->func_aux[0] = (XC(func_type) *)  malloc(sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_LDA_C_PZ, p->nspin);
}


static void 
my_gga_c_p86(const void *p_, const FLOAT *rho, const FLOAT *sigma,
	  FLOAT *e, FLOAT *vrho, FLOAT *vsigma)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  FLOAT dens, zeta, dzdd[2], gdmt, ecunif, vcunif[2];
  FLOAT rs, DD, dDDdzeta, CC, CCinf, dCCdd;
  FLOAT Phi, dPhidd, dPhidgdmt;

  XC(lda_exc_vxc)(p->func_aux[0], 1, rho, &ecunif, vcunif);

  XC(rho2dzeta)(p->nspin, rho, &dens, &zeta);
  if(dens <= 0.0) return;

  dzdd[0] =  (1.0 - zeta)/dens;
  dzdd[1] = -(1.0 + zeta)/dens;
    
  rs = RS(dens);

  /* get gdmt = |nabla n| */
  gdmt = sigma[0];
  if(p->nspin == XC_POLARIZED) gdmt += 2.0*sigma[1] + sigma[2];
  gdmt = SQRT(gdmt);
  if(gdmt < MIN_GRAD) gdmt = MIN_GRAD;


  { /* Equation [1].(4) */ 
    DD       = SQRT(POW(1.0 + zeta, 5.0/3.0) + POW(1.0 - zeta, 5.0/3.0))/M_SQRT2;
    dDDdzeta = 5.0/(3.0*4.0*DD)*(POW(1.0 + zeta, 2.0/3.0) - POW(1.0 - zeta, 2.0/3.0));
  }

  { /* Equation (6) of [1] */
    static const FLOAT alpha = 0.023266, beta = 7.389e-6, gamma = 8.723, delta = 0.472;
    static const FLOAT aa = 0.001667, bb = 0.002568;

    FLOAT rs2 = rs*rs, f1, f2, df1, df2, drsdd;

    f1    = bb + alpha*rs + beta*rs2;
    f2    = 1.0 + gamma*rs + delta*rs2 + 1.0e4*beta*rs*rs2;
    CC    = aa + f1/f2;
    CCinf = aa + bb;

    df1   = alpha + 2.0*beta*rs;
    df2   = gamma + 2.0*delta*rs + 3.0e4*beta*rs2;
    drsdd = -rs/(3.0*dens);
    dCCdd = (df1*f2 - f1*df2)/(f2*f2)*drsdd;
  }

  { /* Equation (9) of [1] */
    static const FLOAT ftilde = 1.745*0.11;

    FLOAT f1, f2, df1, df2;

    f1  = ftilde*(CCinf/CC);
    f2  = POW(dens, -7.0/6.0);
    Phi = f1*gdmt*f2;

    df1 = -f1/(CC)*dCCdd;
    df2 = -7.0/6.0*POW(dens, -13.0/6.0);
    dPhidd    = gdmt*(df1*f2 + f1*df2);
    dPhidgdmt = f1*f2;
  }

  { /* Equation [1].(8) */
    FLOAT gdmt2;
    FLOAT f1, f2, f3, df1, df1dgdmt, df2, df3, df3dgdmt;

    gdmt2 = gdmt*gdmt;

    f1 = exp(-Phi);
    f2 = POW(dens, -4.0/3.0);
    f3 = f1*CC*gdmt2*f2;

    df1      = -f1*dPhidd;
    df1dgdmt = -f1*dPhidgdmt;
    df2      = -4.0/3.0*POW(dens, -7.0/3.0);
    df3      = gdmt2*(df1*CC*f2 + f1*dCCdd*f2 + f1*CC*df2);
    df3dgdmt = CC*f2*(df1dgdmt*gdmt2 + f1*2.0*gdmt);

    *e = ecunif + f3/(DD*dens);

    if(vrho != NULL){
      vrho[0]   = vcunif[0] + (df3 - (f3/DD)*dDDdzeta*dzdd[0])/DD;
      vsigma[0] = df3dgdmt/(DD*2.0*gdmt);

      if(p->nspin == XC_POLARIZED){
	vrho[1]   = vcunif[1] + (df3 - (f3/DD)*dDDdzeta*dzdd[1])/DD;
	vsigma[1] = 2.0*vsigma[0];
	vsigma[2] =     vsigma[0];
      }
    }
  }
}

/* Warning: this is a workaround to support blocks while waiting for the next interface */
static void 
gga_c_p86(const void *p_, int np, const FLOAT *rho, const FLOAT *sigma,
	  FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	  FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  int ip;
  const XC(gga_type) *p = p_;

  for(ip=0; ip<np; ip++){
    my_gga_c_p86(p_, rho, sigma, zk, vrho, vsigma);

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


const XC(func_info_type) XC(func_info_gga_c_p86) = {
  XC_GGA_C_P86,
  XC_CORRELATION,
  "Perdew 86",
  XC_FAMILY_GGA,
  "JP Perdew, Phys. Rev. B 33, 8822 (1986)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  gga_c_p86_init,
  NULL,
  NULL,
  gga_c_p86
};
