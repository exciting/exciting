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

void 
XC(perdew_params)(const XC(gga_type) *gga_p, const FLOAT *rho, const FLOAT *sigma, int order, XC(perdew_t) *pt)
{
  pt->nspin = gga_p->nspin;
  XC(rho2dzeta)(pt->nspin, rho, &(pt->dens), &(pt->zeta));

  if(pt->dens < MIN_DENS) return;

  switch (order){
  case 0:
    XC(lda_exc) (gga_p->func_aux[0], 1, rho, &(pt->ecunif));
    break;
  case 1:
    XC(lda_exc_vxc)(gga_p->func_aux[0], 1, rho, &(pt->ecunif), pt->vcunif);
    break;
  case 2:
    XC(lda)(gga_p->func_aux[0], 1, rho, &(pt->ecunif), pt->vcunif, pt->fcunif, NULL);
    break;
  }

  pt->rs = RS(pt->dens);
  pt->kf = CBRT(3.0*M_PI*M_PI*pt->dens);
  pt->ks = SQRT(4.0*pt->kf/M_PI);

  /* phi is bounded between 2^(-1/3) and 1 */
  pt->phi  = 0.5*(POW(1.0 + pt->zeta, 2.0/3.0) + POW(1.0 - pt->zeta, 2.0/3.0));

  /* get gdmt = |nabla n| */
  pt->gdmt = sigma[0];
  if(pt->nspin == XC_POLARIZED) pt->gdmt += 2.0*sigma[1] + sigma[2];
  pt->gdmt = SQRT(max(pt->gdmt, 0.0));

  pt->t = pt->gdmt/(2.0 * pt->phi * pt->ks * pt->dens);

  if(order > 0)
    pt->drs = pt->dkf = pt->dks = pt->dphi = pt->dt = pt->decunif = 0.0;

  if(order > 1){
    pt->d2rs2 = pt->d2rskf = pt->d2rsks = pt->d2rsphi = pt->d2rst  = pt->d2rsecunif  = 0.0;
                pt->d2kf2  = pt->d2kfks = pt->d2kfphi = pt->d2kft  = pt->d2kfecunif  = 0.0;
		             pt->d2ks2  = pt->d2ksphi = pt->d2kst  = pt->d2ksecunif  = 0.0;
                                          pt->d2phi2  = pt->d2phit = pt->d2phiecunif = 0.0;
                                                        pt->d2t2   = pt->d2tecunif   = 0.0;
                                                                     pt->d2ecunif2   = 0.0;
  }
}

void 
XC(perdew_potentials)(XC(perdew_t) *pt, const FLOAT *rho, FLOAT e_gga, int order,
		      FLOAT *vrho, FLOAT *vsigma,
		      FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  /* alpha = {0->rs, 1->kf, 2->ks, 3->phi, 4->t, 5->ec */
  FLOAT   dalphadd[6][2],   dFdalpha[6];
  FLOAT d2alphadd2[6][3], d2Fdalpha2[6][6];

  FLOAT dzdd[2], dpdz, d2zdd2[3], d2pdz2;
  FLOAT dtdsig, d2tdsig2;
  int is, js, ks, ns;
 
  if(order < 1) return;

  if(pt->nspin == XC_POLARIZED){
    FLOAT aux;

    dpdz    = 0.0;
    /* This is written like this to workaround a problem with xlf compilers */
    if(fabs(1.0 + pt->zeta) >= MIN_DENS){
      aux   = CBRT(1.0 + pt->zeta);
      dpdz += 1.0/(3.0*aux);
    }
    if(fabs(1.0 - pt->zeta) >= MIN_DENS){
      aux   = CBRT(1.0 - pt->zeta);
      dpdz -= 1.0/(3.0*aux);
    }

    dzdd[0] =  (1.0 - pt->zeta)/pt->dens;
    dzdd[1] = -(1.0 + pt->zeta)/pt->dens;
  }else{
    dpdz    = 0.0;
    dzdd[0] = 0.0;
  }

  dFdalpha[0] = pt->drs;
  dFdalpha[1] = pt->dkf;
  dFdalpha[2] = pt->dks;
  dFdalpha[3] = pt->dphi;
  dFdalpha[4] = pt->dt;
  dFdalpha[5] = pt->decunif;

  for(is=0; is<pt->nspin; is++){
    dalphadd[0][is] = -pt->rs/(3.0*pt->dens);
    dalphadd[1][is] =  pt->kf/(3.0*pt->dens);
    dalphadd[2][is] =  pt->ks*dalphadd[1][is]/(2.0*pt->kf);
    dalphadd[3][is] =  dpdz*dzdd[is];
    dalphadd[4][is] = -pt->t*(1.0/pt->dens + dalphadd[2][is]/pt->ks + dalphadd[3][is]/pt->phi);;
    dalphadd[5][is] = (pt->vcunif[is] - pt->ecunif)/pt->dens;
  }

  /* calculate vrho */
  if(vrho != NULL)
    for(is=0; is<pt->nspin; is++){
      if(rho[is] > MIN_DENS){
	int k;
	
	vrho[is] = e_gga;
	for(k=0; k<6; k++)
	  vrho[is] += pt->dens * dFdalpha[k]*dalphadd[k][is];
      }else{
      vrho[is] = 0.0;
      }
    }

  if(pt->gdmt > MIN_GRAD){
    dtdsig  = pt->t/(2.0*pt->gdmt*pt->gdmt);

    if(vrho != NULL){ /* calculate now vsigma */
      vsigma[0] = pt->dens*pt->dt*dtdsig;
      if(pt->nspin == XC_POLARIZED){
	vsigma[1] = 2.0*vsigma[0];
	vsigma[2] =     vsigma[0];
      }
    }
  }

  if(order < 2) return;

  /* first let us sort d2Fdalpha2 in a matrix format */
  d2Fdalpha2[0][0] = pt->d2rs2;
  d2Fdalpha2[0][1] = pt->d2rskf;
  d2Fdalpha2[0][2] = pt->d2rsks;
  d2Fdalpha2[0][3] = pt->d2rst;
  d2Fdalpha2[0][4] = pt->d2rsphi;
  d2Fdalpha2[0][5] = pt->d2rsecunif;

  d2Fdalpha2[1][0] = d2Fdalpha2[0][1];
  d2Fdalpha2[1][1] = pt->d2kf2;
  d2Fdalpha2[1][2] = pt->d2kfks;
  d2Fdalpha2[1][3] = pt->d2kft;
  d2Fdalpha2[1][4] = pt->d2kfphi;
  d2Fdalpha2[1][5] = pt->d2kfecunif;

  d2Fdalpha2[2][0] = d2Fdalpha2[0][2];
  d2Fdalpha2[2][1] = d2Fdalpha2[1][2];
  d2Fdalpha2[2][2] = pt->d2ks2;
  d2Fdalpha2[2][3] = pt->d2kst;
  d2Fdalpha2[2][4] = pt->d2ksphi;
  d2Fdalpha2[2][5] = pt->d2ksecunif;

  d2Fdalpha2[3][0] = d2Fdalpha2[0][3];
  d2Fdalpha2[3][1] = d2Fdalpha2[1][3];
  d2Fdalpha2[3][2] = d2Fdalpha2[2][3];
  d2Fdalpha2[3][3] = pt->d2phi2;
  d2Fdalpha2[3][4] = pt->d2phit;
  d2Fdalpha2[3][5] = pt->d2phiecunif;

  d2Fdalpha2[4][0] = d2Fdalpha2[0][4];
  d2Fdalpha2[4][1] = d2Fdalpha2[1][4];
  d2Fdalpha2[4][2] = d2Fdalpha2[2][4];
  d2Fdalpha2[4][3] = d2Fdalpha2[3][4];
  d2Fdalpha2[4][4] = pt->d2t2;
  d2Fdalpha2[4][5] = pt->d2tecunif;

  d2Fdalpha2[5][0] = d2Fdalpha2[0][5];
  d2Fdalpha2[5][1] = d2Fdalpha2[1][5];
  d2Fdalpha2[5][2] = d2Fdalpha2[2][5];
  d2Fdalpha2[5][3] = d2Fdalpha2[3][5];
  d2Fdalpha2[5][4] = d2Fdalpha2[4][5];
  d2Fdalpha2[5][5] = pt->d2ecunif2;

  /* now we sort d2alphadd2 */
  if(pt->nspin == XC_POLARIZED){
    FLOAT aux;

    d2pdz2 = 0.0;
    if(fabs(1.0 + pt->zeta) >= MIN_DENS){
      aux     = POW(1.0 + pt->zeta, 4.0/3.0);
      d2pdz2 += -1.0/(9.0*aux);
    }
    if(fabs(1.0 - pt->zeta) >= MIN_DENS){
      aux     = POW(1.0 - pt->zeta, 4.0/3.0);
      d2pdz2 += -1.0/(9.0*aux);
    }

    d2zdd2[0] = -2.0*dzdd[0]/pt->dens;
    d2zdd2[1] =  2.0*pt->zeta/(pt->dens*pt->dens);
    d2zdd2[2] = -2.0*dzdd[1]/pt->dens;
  }else{
    d2pdz2    = 0.0;
    d2zdd2[0] = 0.0;
  }

  ns = (pt->nspin == XC_UNPOLARIZED) ? 0 : 2;
  for(ks=0; ks<=ns; ks++){
    is = (ks == 0 || ks == 1) ? 0 : 1;
    js = (ks == 0           ) ? 0 : 1;

    d2alphadd2[0][ks] =  4.0/9.0*pt->rs/(pt->dens*pt->dens);

    d2alphadd2[1][ks] = -2.0/9.0*pt->kf/(pt->dens*pt->dens);

    d2alphadd2[2][ks] =  pt->ks/(2.0*pt->kf)*
      (d2alphadd2[1][ks] - dalphadd[1][is]*dalphadd[1][js]/(2.0*pt->kf));

    d2alphadd2[3][ks] =  d2pdz2*dzdd[is]*dzdd[js] + dpdz*d2zdd2[ks];

    d2alphadd2[4][ks] =  pt->t *
      (+2.0/(pt->dens*pt->dens)
       +2.0/(pt->ks*pt->ks)   *(dalphadd[2][is] * dalphadd[2][js])
       +2.0/(pt->phi*pt->phi) *(dalphadd[3][is] * dalphadd[3][js])
       +1.0/(pt->dens*pt->ks) *(dalphadd[2][is] + dalphadd[2][js])
       +1.0/(pt->dens*pt->phi)*(dalphadd[3][is] + dalphadd[3][js])
       +1.0/(pt->ks*pt->phi)  *(dalphadd[2][is]*dalphadd[3][js] + dalphadd[2][js]*dalphadd[3][is])
       -1.0/(pt->ks)*d2alphadd2[2][ks] -1.0/(pt->phi)*d2alphadd2[3][ks]);

    d2alphadd2[5][ks] = pt->fcunif[ks]/pt->dens -
      (pt->vcunif[is] + pt->vcunif[js] - 2.0*pt->ecunif)/(pt->dens*pt->dens);
  }

  for(ks=0; ks<=ns; ks++){
    int j, k;

    is = (ks == 0 || ks == 1) ? 0 : 1;
    js = (ks == 0           ) ? 0 : 1;

    v2rho2[ks] = 0.0;

    for(j=0; j<6; j++){
      v2rho2[ks] += dFdalpha[j]*(dalphadd[j][is] + dalphadd[j][js]);
      v2rho2[ks] += pt->dens * dFdalpha[j]*d2alphadd2[j][ks];

      for(k=0; k<6; k++)
	v2rho2[ks] +=  pt->dens * d2Fdalpha2[j][k]*dalphadd[j][is]*dalphadd[k][js];
    }
  }

  /* now we handle v2rhosigma */
  if(pt->gdmt > MIN_GRAD){
    for(is=0; is<pt->nspin; is++){
      int j;
      ks = (is == 0) ? 0 : 5;
      
      v2rhosigma[ks] = dFdalpha[4]*dtdsig;

      for(j=0; j<6; j++)
	v2rhosigma[ks] += pt->dens * d2Fdalpha2[4][j]*dalphadd[j][is]*dtdsig;

      v2rhosigma[ks] += pt->dens * dFdalpha[4]*dalphadd[4][is]/(2.0*pt->gdmt*pt->gdmt);
    }

    if(pt->nspin == XC_POLARIZED){
      v2rhosigma[1] = 2.0*v2rhosigma[0];
      v2rhosigma[2] =     v2rhosigma[0];
      v2rhosigma[3] =     v2rhosigma[5];
      v2rhosigma[4] = 2.0*v2rhosigma[5];
    }

    /* now wwe take care of v2sigma2 */
    d2tdsig2 = -dtdsig/(2.0*pt->gdmt*pt->gdmt);
    v2sigma2[0] = pt->dens*(pt->d2t2*dtdsig*dtdsig + pt->dt*d2tdsig2);
    if(pt->nspin == XC_POLARIZED){
      v2sigma2[1] = 2.0*v2sigma2[0]; /* aa_ab */
      v2sigma2[2] =     v2sigma2[0]; /* aa_bb */
      v2sigma2[3] = 4.0*v2sigma2[0]; /* ab_ab */
      v2sigma2[4] = 2.0*v2sigma2[0]; /* ab_bb */
      v2sigma2[5] =     v2sigma2[0]; /* bb_bb */
    }
  }
  
}
