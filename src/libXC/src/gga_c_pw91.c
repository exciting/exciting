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
#include <math.h>
#include "util.h"

#define XC_GGA_C_PW91 134 /* Perdew & Wang 91 */

/* My implementation differs from the repository due to a constant in the
   PW92 part of the functional. However, if the line (code from repository)
   
   t55 = 0.3799575D-1*t33*t44*t52

   is replaced by

   t55 = 0.3799574853701528D-1*t33*t44*t52

   both results agree. I already mailed Huub van Dam to try to clarify the problem.
*/
static FLOAT pw91_nu, pw91_beta;
static const FLOAT
  pw91_C_c0  = 4.235e-3, 
  pw91_alpha = 0.09;

static void gga_c_pw91_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_LDA_C_PW, p->nspin);

  pw91_nu   = 16.0/M_PI * CBRT(3.0*M_PI*M_PI);
  pw91_beta = pw91_nu*pw91_C_c0;
}


static void gga_c_pw91_end(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  XC(func_end)(p->func_aux[0]);
  free(p->func_aux[0]);
  free(p->func_aux);
}


static void
A_eq14(FLOAT ec, FLOAT g, FLOAT *A, FLOAT *dec, FLOAT *dg)
{
  FLOAT g2, g3, dd;

  g2 = g*g;
  g3 = g*g2;

  dd = -2.0*pw91_alpha*ec/(g3*pw91_beta*pw91_beta);
  dd = exp(dd);

  *A   = (2.0*pw91_alpha/pw91_beta) / (dd - 1.0);

  *dec = -(*A)/(dd - 1.0) * dd * (-2.0*pw91_alpha/(g3*pw91_beta*pw91_beta));
  *dg  = -(*A)/(dd - 1.0) * dd * (+6.0*pw91_alpha*ec/(g*g3*pw91_beta*pw91_beta));
}

static void
H0_eq13(FLOAT   ec, FLOAT   g, FLOAT   t, FLOAT *H0,
	FLOAT *dec, FLOAT *dg, FLOAT *dt)
{
  FLOAT A, dAdec, dAdg;
  FLOAT g3, t2, t4, n0, d0, dd, dA;

  A_eq14(ec, g, &A, &dAdec, &dAdg);

  g3 = g*g*g;
  t2 = t*t;
  t4 = t2*t2;

  n0 = t2 + A*t4;
  d0 = 1.0 + A*t2 + A*A*t4;

  *H0 = g3 * pw91_beta*pw91_beta/(2.0*pw91_alpha) *
    log(1.0 + 2.0*pw91_alpha/pw91_beta * n0/d0);

  dd = d0*(pw91_beta + n0*(2.0*pw91_alpha + A*pw91_beta));
  dA = -A*g3*t2*t4*(2.0 + A*t2)*pw91_beta*pw91_beta/dd;

  *dec = dA * dAdec;
  *dg  = (*H0)*3.0/g + dA*dAdg;
  *dt  = 2.0*g3*t*(1.0 + 2.0*A*t2)*pw91_beta*pw91_beta/dd;
}


/* pade parametrized form of C-xc found in
   M Rasolt & DJW Geldart, Phys. Rev. B 34, 1325 (1986)
*/
static inline void 
Rasold_Geldart_C_xc(FLOAT rs, FLOAT *C_xc, FLOAT *drs)
{
  const FLOAT 
    a[3] = {2.568, 23.266, 0.007389},
    b[3] = {1.0, 8.723, 0.472};
  
  FLOAT d0, d1, n0, n1;

  n0 = (a[0] + rs*(a[1] + rs*a[2]));
  d0 =  b[0] + rs*(b[1] + rs*(b[2] + 10.0*rs*a[2]));

  n1 = a[1] + 2.0*rs*a[2];
  d1 = b[1] + 2.0*rs*b[2] + 10.0*3.0*rs*rs*a[2];

  *C_xc = n0/(1000.0*d0);
  *drs  = (n1*d0 - n0*d1)/(1000.0*d0*d0);
}


static void 
H1_eq15(FLOAT   rs, FLOAT   g, FLOAT   t, FLOAT   ks, FLOAT   kf, FLOAT *H1,
	FLOAT *drs, FLOAT *dg, FLOAT *dt, FLOAT *dks, FLOAT *dkf)
{
  const FLOAT C_xc0 = 2.568e-3, C_x = -0.001667;

  FLOAT g3, g4, t2, kf2, ks2, dd1, dd2;
  FLOAT C_xc, dC_xc;

  g3  = g*g*g;
  g4  = g3*g;
  t2  = t*t;
  ks2 = ks*ks;
  kf2 = kf*kf;

  dd1 = -100.0 * g4 * (ks2/kf2) * t2;
  dd1 = exp(dd1);

  Rasold_Geldart_C_xc(rs, &C_xc, &dC_xc);
  dd2 = C_xc - C_xc0 - 3.0*C_x/7.0;

  *H1  = pw91_nu * dd2 * g3 * t2 * dd1;

  *drs = pw91_nu * dC_xc * g3 * t2 * dd1;
  *dg  = (*H1) * (3.0/g - 100.0* 4.0*g3 *(ks2/kf2)*t2); /* g can not be zero */
  *dt  = pw91_nu * dd2 * g3 * dd1 * 
    (2.0*t - t2*100.0*g4*(ks2/kf2)* 2.0*t);
  *dks = (*H1) * (-100.0*g4*( 2.0*ks /kf2)*t2);
  *dkf = (*H1) * (+100.0*g4*( 2.0*ks2/(kf*kf2) )*t2);
}


static inline void 
ec_eq9(FLOAT   ec, FLOAT   rs, FLOAT   t, FLOAT   g, FLOAT   ks, FLOAT   kf, FLOAT  *ec_gga,
       FLOAT *dec, FLOAT *drs, FLOAT *dt, FLOAT *dg, FLOAT *dks, FLOAT *dkf)
{
  FLOAT H0, dH0dec, dH0dg, dH0dt;
  FLOAT H1, dH1drs, dH1dg, dH1dt, dH1dks, dH1dkf;

  H0_eq13(ec, g, t, &H0, 
	  &dH0dec, &dH0dg, &dH0dt);
  H1_eq15(rs, g, t, ks, kf, &H1,
	  &dH1drs, &dH1dg, &dH1dt, &dH1dks, &dH1dkf);

  *ec_gga = ec + H0 + H1;
  *dec    = 1.0 + dH0dec;
  *drs    = dH1drs;
  *dt     = dH0dt + dH1dt;
  *dg     = dH0dg + dH1dg;
  *dks    = dH1dks;
  *dkf    = dH1dkf;
}

static void 
my_gga_c_pw91(const void *p_, const FLOAT *rho, const FLOAT *sigma,
	   FLOAT *e, FLOAT *vrho, FLOAT *vsigma,
	   FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  int order;
  XC(gga_type) *p = (XC(gga_type) *)p_;
  XC(perdew_t) pt;

  order = 0;
  if(vrho   != NULL) order = 1;

  XC(perdew_params)(p, rho, sigma, order, &pt);
  if(pt.dens < MIN_DENS) return;

  ec_eq9(pt.ecunif, pt.rs, pt.t, pt.phi, pt.ks, pt.kf, e,
	 &pt.decunif, &pt.drs, &pt.dt, &pt.dphi, &pt.dks, &pt.dkf);

  XC(perdew_potentials)(&pt, rho, *e, order, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);
}

/* Warning: this is a workaround to support blocks while waiting for the next interface */
static void 
gga_c_pw91(const void *p_, int np, const FLOAT *rho, const FLOAT *sigma,
	  FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	  FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  int ip;
  const XC(gga_type) *p = p_;

  for(ip=0; ip<np; ip++){
    my_gga_c_pw91(p_, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);

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

const XC(func_info_type) XC(func_info_gga_c_pw91) = {
  XC_GGA_C_PW91,
  XC_CORRELATION,
  "Perdew & Wang 91",
  XC_FAMILY_GGA,
  "JP Perdew, JA Chevary, SH Vosko, KA Jackson, MR Pederson, DJ Singh, and C Fiolhais, Phys. Rev. B 46, 6671 (1992)\n"
  "JP Perdew, JA Chevary, SH Vosko, KA Jackson, MR Pederson, DJ Singh, and C Fiolhais, Phys. Rev. B 48, 4978(E) (1993)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  gga_c_pw91_init,
  gga_c_pw91_end,
  NULL,            /* this is not an LDA                   */
  gga_c_pw91,
};
