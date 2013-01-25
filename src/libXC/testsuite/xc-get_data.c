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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <xc.h>

typedef struct {
  int functional;
  int nspin;

  /* Input */
  double rho[2];         /* rhoa, rhob */
  double sigma[3];       /* sigmaaa, sigmaab, sigmabb */
  double lapl[2];        /* lapla, laplb */
  double tau[2];         /* taua, taub */

  /* Energy */
  double zk;             /* energy density per unit particle */

  /* First derivatives */
  double vrho[2];        /* vrhoa, vrhob */
  double vsigma[3];      /* vsigmaaa, vsigmaab, vsigmabb */
  double vlapl[2];       /* vlapla, vlaplb */
  double vtau[2];        /* vtaua, vtaub */

  /* Second derivatives */
  double v2rho2[3];      /* v2rhoa2, v2rhoab, v2rhob2 */
  double v2rhosigma[6];  /* v2rhoasigmaaa, v2rhoasigmaab, v2rhoasigmabb
			    v2rhobsigmaaa, v2rhobsigmaab, v2rhobsigmabb */
  double v2rholapl[3];   /* */
  double v2rhotau[3];    /* */
  double v2sigma2[6];    /* v2sigmaaa2, v2sigmaaaab, v2sigmaaabb
			    v2sigmaab2, v2sigmaabbb, v2sigmabb2 */
  double v2sigmalapl[6]; /* v2sigmaaalapla, v2sigmaaalaplb,
			    v2sigmaablapla, v2sigmaablaplb,
			    v2sigmabblapla, v2sigmabblaplb */
  double v2sigmatau[6];  /* v2sigmaaataua, v2sigmaaataub, 
			    v2sigmaabtaua, v2sigmaabtaub,
			    v2sigmabbtaua, v2sigmabbtaub */
  double v2lapl2[3];     /* v2lapla2, v2laplab, v2laplb2 */
  double v2lapltau[3];   /* */
  double v2tau2[3];      /* v2taua2, v2tauab, v2taub2 */

  /* Third derivatives */
  double v3rho3[4];      /* v3rhoaaa, v3rhoaab, v3rhoabb, v3rhobbb */

} xc_values_type;

/*----------------------------------------------------------*/
void init_values(xc_values_type *xc_values, char *argv[])
{
  int i;

  xc_values->functional = atoi(argv[1]);
  xc_values->nspin      = atoi(argv[2]);
  xc_values->rho[0]     = atof(argv[3]);
  xc_values->rho[1]     = atof(argv[4]);
  xc_values->sigma[0]   = atof(argv[5]);
  xc_values->sigma[1]   = atof(argv[6]);
  xc_values->sigma[2]   = atof(argv[7]);
  xc_values->lapl[0]    = atof(argv[8]);
  xc_values->lapl[1]    = atof(argv[9]);
  xc_values->tau[0]     = atof(argv[10]);
  xc_values->tau[1]     = atof(argv[11]);

  xc_values->zk = 0;

  for(i=0; i<2; i++){
    xc_values->vrho[i]  = 0;
    xc_values->vlapl[i] = 0;
    xc_values->vtau[i]  = 0;
  }

  for(i=0; i<3; i++){
    xc_values->vsigma[i]    = 0;
    xc_values->v2rho2[i]    = 0;
    xc_values->v2lapl2[i]   = 0;
    xc_values->v2tau2[i]    = 0;
    xc_values->v2rholapl[i] = 0;
    xc_values->v2rhotau[i]  = 0;
    xc_values->v2lapltau[i] = 0;
  }

  for(i=0; i<4; i++){
    xc_values->v3rho3[i]  = 0;
  }

  for(i=0; i<6; i++){
    xc_values->v2rhosigma[i]  = 0;
    xc_values->v2sigma2[i]    = 0;
    xc_values->v2sigmalapl[i] = 0;
    xc_values->v2sigmatau[i]  = 0;
  }
}


/*----------------------------------------------------------*/
void print_values(xc_values_type *xc)
{
  //int family = xc_family_from_id(xc->functional, NULL, NULL);

  printf(" rhoa= %#0.2E rhob= %#0.2E sigmaaa= %#0.2E sigmaab= %#0.2E sigmabb= %#0.2E lapla= %#0.2E laplb= %#0.2E taua= %#0.2E taub= %#0.2E\n\n",
	 xc->rho[0], xc->rho[1],
	 xc->sigma[0], xc->sigma[1], xc->sigma[2],
	 xc->lapl[0], xc->lapl[1],
	 xc->tau[0], xc->tau[1]);
  printf(" zk            = %#19.12E\n\n",
	 xc->zk);
  printf(" vrhoa         = %#19.12E\n"
	 " vrhob         = %#19.12E\n"
	 " vsigmaaa      = %#19.12E\n"
	 " vsigmaab      = %#19.12E\n"
	 " vsigmabb      = %#19.12E\n"
	 " vlapla        = %#19.12E\n"
	 " vlaplb        = %#19.12E\n"
	 " vtaua         = %#19.12E\n"
	 " vtaub         = %#19.12E\n\n",
	 xc->vrho[0], xc->vrho[1],
	 xc->vsigma[0], xc->vsigma[1], xc->vsigma[2],
	 xc->vlapl[0], xc->vlapl[1],
	 xc->vtau[0], xc->vtau[1]);

  printf(" v2rhoa2       = %#19.12E\n"
	 " v2rhoab       = %#19.12E\n"
	 " v2rhob2       = %#19.12E\n"
	 " v2rhoasigmaaa = %#19.12E\n"
	 " v2rhoasigmaab = %#19.12E\n"
	 " v2rhoasigmabb = %#19.12E\n"
	 " v2rhobsigmaaa = %#19.12E\n"
	 " v2rhobsigmaab = %#19.12E\n"
	 " v2rhobsigmabb = %#19.12E\n"
	 " v2sigmaaa2    = %#19.12E\n"
	 " v2sigmaaaab   = %#19.12E\n"
	 " v2sigmaaabb   = %#19.12E\n"
	 " v2sigmaab2    = %#19.12E\n"
	 " v2sigmaabbb   = %#19.12E\n"
	 " v2sigmabb2    = %#19.12E\n\n",
	 xc->v2rho2[0], xc->v2rho2[1], xc->v2rho2[2],
	 xc->v2rhosigma[0], xc->v2rhosigma[1], xc->v2rhosigma[2],
	 xc->v2rhosigma[3], xc->v2rhosigma[4], xc->v2rhosigma[5],
	 xc->v2sigma2[0], xc->v2sigma2[1], xc->v2sigma2[2],
	 xc->v2sigma2[3], xc->v2sigma2[4], xc->v2sigma2[5]
	 );
  printf(" v3rhoa3  = %#19.12E\n"
	 " v2rhoaab = %#19.12E\n"
	 " v2rhoabb = %#19.12E\n"
	 " v2rhob3  = %#19.12E\n\n",
	 xc->v3rho3[0], xc->v3rho3[1], xc->v3rho3[2], xc->v3rho3[3]
	 );

}


/*----------------------------------------------------------*/
int main(int argc, char *argv[])
{
  xc_values_type xc;
  xc_func_type func;
  const xc_func_info_type *info;

  FLOAT *pzk          = NULL;
  FLOAT *pvrho        = NULL;
  FLOAT *pvsigma      = NULL;
  FLOAT *pvlapl       = NULL;
  FLOAT *pvtau        = NULL;
  FLOAT *pv2rho2      = NULL;
  FLOAT *pv2rhosigma  = NULL;
  FLOAT *pv2rholapl   = NULL;
  FLOAT *pv2rhotau    = NULL;
  FLOAT *pv2sigma2    = NULL;
  FLOAT *pv2sigmalapl = NULL;
  FLOAT *pv2sigmatau  = NULL;
  FLOAT *pv2lapl2     = NULL;
  FLOAT *pv2lapltau   = NULL;
  FLOAT *pv2tau2      = NULL;
  FLOAT *pv3rho3      = NULL;

  if(argc != 12){
    printf("Usage:\n%s funct pol rhoa rhob sigmaaa sigmaab sigmabb lapla laplb taua taub\n", argv[0]);
    return 1;
  }

  init_values(&xc, argv);

  if(xc.nspin == 1){
    xc.rho[0]   += xc.rho[1];
    xc.sigma[0] += 2.0*xc.sigma[1] + xc.sigma[2];
    xc.lapl[0]  += xc.lapl[1];
    xc.tau[0]   += xc.tau[1];
  }

  if(xc_func_init(&func, xc.functional, xc.nspin) != 0){
    fprintf(stderr, "Functional '%d' not found\n", xc.functional);
    exit(1);  
  }
  info = func.info;

  if(info->flags & XC_FLAGS_HAVE_EXC){
    pzk = &xc.zk;
  }
  if(info->flags & XC_FLAGS_HAVE_VXC){
    pvrho   = xc.vrho;
    pvsigma = xc.vsigma;
    pvlapl  = xc.vlapl;
    pvtau   = xc.vtau;
  }
  if(info->flags & XC_FLAGS_HAVE_FXC){
    pv2rho2      = xc.v2rho2;
    pv2rhosigma  = xc.v2rhosigma;
    pv2rholapl   = xc.v2rholapl;
    pv2rhotau    = xc.v2rhotau;
    pv2sigma2    = xc.v2sigma2;
    pv2sigmalapl = xc.v2sigmalapl;
    pv2sigmatau  = xc.v2sigmatau;
    pv2lapl2     = xc.v2lapl2;
    pv2lapltau   = xc.v2lapltau;
    pv2tau2      = xc.v2tau2;
  }
  if(info->flags & XC_FLAGS_HAVE_KXC){
    pv3rho3 = xc.v3rho3;
  }

  switch(func.info->family)
    {
    case XC_FAMILY_LDA:
      xc_lda(&func, 1, xc.rho, pzk, pvrho, pv2rho2, pv3rho3);
      break;
    case XC_FAMILY_GGA:
    case XC_FAMILY_HYB_GGA:
      xc_gga(&func, 1, xc.rho, xc.sigma,
	     pzk, pvrho, pvsigma, pv2rho2, pv2rhosigma, pv2sigma2);
      break;
    case XC_FAMILY_MGGA:
      xc_mgga(&func, 1, xc.rho, xc.sigma, xc.lapl, xc.tau,
	      pzk, pvrho, pvsigma, pvlapl, pvtau,
	      pv2rho2, pv2sigma2, pv2lapl2, pv2tau2, pv2rhosigma, pv2rholapl, pv2rhotau, pv2sigmalapl, pv2sigmatau, pv2lapltau);
      break;
    }

  xc_func_end(&func);

  if(xc.nspin == 1){
    xc.zk *= xc.rho[0];
  }else{
    xc.zk *= (xc.rho[0] + xc.rho[1]);
  }


  print_values(&xc);

  return 0;
}

