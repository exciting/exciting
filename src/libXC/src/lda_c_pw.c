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

/************************************************************************
 Correlation energy per-particle and potential of a HEG as parameterized 
 by 
   J.P. Perdew & Y. Wang
   Ortiz & Ballone

Note that the PW modified, corresponds to the version of PW used in the 
original PBE routine. This amounts to adding some more digits in some of
the constants of PW.
************************************************************************/

#define XC_LDA_C_PW     12   /* Perdew & Wang                */
#define XC_LDA_C_PW_MOD 13   /* Perdew & Wang (Modified)     */
#define XC_LDA_C_OB_PW  14   /* Ortiz & Ballone (PW)         */


/* Function g defined by Eq. 10 of the original paper,
   and it's derivative with respect to rs, Eq. A5 */
static void g(int func, int order, int k, FLOAT *rs, 
	      FLOAT *f, FLOAT *dfdrs, FLOAT *d2fdrs2, FLOAT *d3fdrs3)
{
  static FLOAT a[3][3]     = 
    {
      {0.031091,  0.015545,   0.016887},    /* PW */
      {0.0310907, 0.01554535, 0.0168869},   /* PW (modified) */
      {0.031091,  0.015545,   0.016887}     /* OB */
  }; 
  static FLOAT alpha[3][3] = 
    {
      {0.21370,  0.20548,  0.11125},    /* PW */
      {0.21370,  0.20548,  0.11125},    /* PW (modified) */
      {0.026481, 0.022465, 0.11125}     /* OB */
    };
  static FLOAT beta[3][3][4] = {
    {
      { 7.5957,  3.5876,   1.6382,  0.49294}, /* PW */
      {14.1189,  6.1977,   3.3662,  0.62517},
      {10.357,   3.6231,   0.88026, 0.49671}
    },{
      { 7.5957,  3.5876,   1.6382,  0.49294}, /* PW (modified) */
      {14.1189,  6.1977,   3.3662,  0.62517},
      {10.357,   3.6231,   0.88026, 0.49671}
    },{
      { 7.5957,  3.5876,  -0.46647, 0.13354}, /* OB */
      {14.1189,  6.1977,  -0.56043, 0.11313},
      {10.357,   3.6231,   0.88026, 0.49671}
    }};
  
  FLOAT q0, dq0, q1, dq1, q2, d2q1, d3q1, aux1;
  
  q0  = -2.0*a[func][k]*(1.0 + alpha[func][k]*rs[1]);
  q1  =  2.0*a[func][k];
  q1 *= beta[func][k][0]*rs[0] + beta[func][k][1]*rs[1] + 
    beta[func][k][2]*rs[0]*rs[1] + beta[func][k][3]*rs[2];
  q2  = log(1.0 + 1.0/q1);

  /* the function */
  *f = q0*q2;
  
  if(order < 1) return; /* nothing else to do */

  aux1 = q1*(1.0 + q1);

  /* and now the derivative */
  dq0 = -2.0*a[func][k]*alpha[func][k];
  dq1 = a[func][k]*(beta[func][k][0]/rs[0] + 2.0*beta[func][k][1] + 
		    3.0*beta[func][k][2]*rs[0] + 4.0*beta[func][k][3]*rs[1]);

  *dfdrs = dq0*q2 - q0*dq1/aux1;

  if(order < 2) return;

  d2q1 = a[func][k]*(-beta[func][k][0]/(2.0*rs[0]*rs[1]) +
		     3.0*beta[func][k][2]/(2.0*rs[0]) + 4.0*beta[func][k][3]);

  *d2fdrs2 = 1.0/aux1*(-2*dq0*dq1 - q0*d2q1 + q0*(2.0*q1 + 1.0)*dq1*dq1/aux1);

  if(order < 3) return;
  
  d3q1 = (3.0/4.0)*a[func][k]*(beta[func][k][0]/(rs[0]*rs[2]) - beta[func][k][2]/(rs[0]*rs[1]));

  *d3fdrs3  =  2.0*q0*(2.0 + 3.0*q1)*dq1*dq1*dq1;
  *d3fdrs3 += -3.0*aux1*dq1*(dq0*dq1 + q0*d2q1);
  *d3fdrs3 += (1.0 + q1)*(1.0 + q1)*
    (-6.0*q0*dq1*dq1*dq1 + 6.0*q1*dq1*(dq0*dq1 + q0*d2q1)
     -q1*q1*(3.0*dq0*d2q1 + q0*d3q1));
  *d3fdrs3 /= aux1*aux1*aux1;
}


/* the functional */
static inline void 
func(const XC(lda_type) *p, XC(lda_rs_zeta) *r)
{
  static FLOAT fz20[3] = {
    1.709921,                           /* PW */
    1.709920934161365617563962776245,   /* PW (modified) */
    1.709921                            /* OB */
  };

  int func;
  FLOAT ecp, vcp, fcp, kcp;
  FLOAT ecf, vcf, fcf, kcf;
  FLOAT alpha, dalpha, d2alpha, d3alpha;
  FLOAT z2, z3, z4, fz, dfz, d2fz, d3fz;  

  func = p->info->number - XC_LDA_C_PW;
  assert(func==0 || func==1 || func==2);
  
  /* ec(rs, 0) */
  g(func, r->order, 0, r->rs, &ecp, &vcp, &fcp, &kcp);
  
  if(p->nspin == XC_UNPOLARIZED)
    r->zk = ecp;
  else{
    /* get ferromagnetic values */
    g(func, r->order, 1, r->rs, &ecf, &vcf, &fcf, &kcf);

    /* get -alpha_c */
    g(func, r->order, 2, r->rs, &alpha, &dalpha, &d2alpha, &d3alpha);
    alpha *= -1.0;

    fz  = FZETA(r->zeta);
    z2  = r->zeta*r->zeta;
    z3  = r->zeta*z2;
    z4  = r->zeta*z3;
    r->zk = ecp + z4*fz*(ecf - ecp - alpha/fz20[func]) + fz*alpha/fz20[func];
  }

  if(r->order < 1) return;

  if(p->nspin == XC_UNPOLARIZED)
    r->dedrs = vcp;
  else{
    dalpha *= -1.0;

    dfz = DFZETA(r->zeta);
    r->dedrs = vcp + z4*fz*(vcf - vcp - dalpha/fz20[func]) + fz*dalpha/fz20[func];
    r->dedz  = (4.0*z3*fz + z4*dfz)*(ecf - ecp - alpha/fz20[func])
      + dfz*alpha/fz20[func];
  }

  if(r->order < 2) return;

  if(p->nspin == XC_UNPOLARIZED)
    r->d2edrs2 = fcp;
  else{
    d2alpha *= -1.0;

    d2fz = D2FZETA(r->zeta);
    r->d2edrs2 = fcp + z4*fz*(fcf - fcp - d2alpha/fz20[func]) + fz*d2alpha/fz20[func];
    r->d2edrsz = (4.0*z3*fz + z4*dfz)*(vcf - vcp - dalpha/fz20[func])
      + dfz*dalpha/fz20[func];
    r->d2edz2  = (12.0*z2*fz + 8.0*z3*dfz + z4*d2fz)*(ecf - ecp - alpha/fz20[func])
      + d2fz*alpha/fz20[func];
  }

  if(r->order < 3) return;

  if(p->nspin == XC_UNPOLARIZED)
    r->d3edrs3 = kcp;
  else{
    d3alpha *= -1.0;

    d3fz = D3FZETA(r->zeta);

    r->d3edrs3  = kcp + z4*fz*(kcf - kcp - d3alpha/fz20[func]) + fz*d3alpha/fz20[func];
    r->d3edrs2z = (4.0*z3*fz + z4*dfz)*(fcf - fcp - d2alpha/fz20[func])
      + dfz*d2alpha/fz20[func];
    r->d3edrsz2 = (12.0*z2*fz + 8.0*z3*dfz + z4*d2fz)*(vcf - vcp - dalpha/fz20[func])
      + d2fz*dalpha/fz20[func];
    r->d3edz3   = (24.0*r->zeta*fz + 36.0*z2*dfz + 12*z3*d2fz + z4*d3fz)*(ecf - ecp - alpha/fz20[func])
      + d3fz*alpha/fz20[func];
  }
}

#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_c_pw) = {
  XC_LDA_C_PW,
  XC_CORRELATION,
  "Perdew & Wang",
  XC_FAMILY_LDA,
  "JP Perdew and Y Wang, Phys. Rev. B 45, 13244 (1992)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  NULL,     /* init */
  NULL,     /* end  */
  work_lda, /* lda  */
};

const XC(func_info_type) XC(func_info_lda_c_pw_mod) = {
  XC_LDA_C_PW_MOD,
  XC_CORRELATION,
  "Perdew & Wang (modified)",
  XC_FAMILY_LDA,
  "JP Perdew and Y Wang, Phys. Rev. B 45, 13244 (1992)\n"
  "Added extra digits to some constants as in the PBE routine\n"
  "http://dft.rutgers.edu/pubs/PBE.asc",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  NULL,     /* init */
  NULL,     /* end  */
  work_lda, /* lda  */
};

const XC(func_info_type) XC(func_info_lda_c_ob_pw) = {
  XC_LDA_C_OB_PW,
  XC_CORRELATION,
  "Ortiz & Ballone (PW parametrization)",
  XC_FAMILY_LDA,
  "G Ortiz and P Ballone, Phys. Rev. B 50, 1391 (1994)\n"
  "G Ortiz and P Ballone, Phys. Rev. B 56, 9970(E) (1997)\n"
  "JP Perdew and Y Wang, Phys. Rev. B 45, 13244 (1992)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  NULL,     /* init */
  NULL,     /* end  */
  work_lda, /* lda  */
};
