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

#define XC_GGA_C_LYP  131 /* Lee, Yang & Parr */

typedef struct{
  FLOAT A, B, c, d;
} gga_c_lyp_params;


void XC(gga_c_lyp_init)(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;
  gga_c_lyp_params *params;

  assert(p->params == NULL);

  p->params = malloc(sizeof(gga_c_lyp_params));
  params = (gga_c_lyp_params *) (p->params);

  /* values of constants in standard LYP functional */
  XC(gga_c_lyp_set_params_)(p, 0.04918, 0.132, 0.2533, 0.349);
}


void XC(gga_c_lyp_set_params)(XC(func_type) *p, FLOAT A, FLOAT B, FLOAT c, FLOAT d)
{
  assert(p != NULL && p->gga != NULL);
  XC(gga_c_lyp_set_params_)(p->gga, A, B, c, d);
}


void XC(gga_c_lyp_set_params_)(XC(gga_type) *p, FLOAT A, FLOAT B, FLOAT c, FLOAT d)
{
  gga_c_lyp_params *params;

  assert(p->params != NULL);
  params = (gga_c_lyp_params *) (p->params);

  params->A = A;
  params->B = B;
  params->c = c;
  params->d = d;
}


static void 
my_gga_c_lyp(const void *p_, const FLOAT *rho_, const FLOAT *sigma_,
	     FLOAT *e, FLOAT *vrho, FLOAT *vsigma)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;
  gga_c_lyp_params *params;

  static FLOAT ee = 36.462398978764767321; /* ee = 8*2^(2/3)*e */

  FLOAT rho[2], sigma[2], rhot, sigmat;
  FLOAT AA, BB, cc, dd; /* sortcuts for parameters */
  FLOAT sfact, rhot13, rhot43, rho83[2], ZZ, delta, omega;
  FLOAT dZZdr, ddeltadr, domegadr;
  FLOAT t1, t2, t3, t4, t5, t6;
  int is;

  assert(p->params != NULL);
  params = (gga_c_lyp_params *)(p->params);

  AA = params->A;
  BB = params->B;
  cc = params->c;
  dd = params->d;

  /* convert input to spin-polarized case */
  if(p->nspin == XC_POLARIZED){
    sfact    = 1.0;
    rho[0]   = rho_[0];
    rho[1]   = rho_[1];
    sigma[0] = sigma_[0];
    sigma[1] = sigma_[2];
    rhot     = rho_[0] + rho_[1];
    sigmat   = sigma_[0] + 2.0*sigma_[1] + sigma_[2];
  }else{
    sfact    = 2.0;
    rho[0]   = rho_[0]/2.0;
    rho[1]   = rho_[0]/2.0;
    sigma[0] = sigma_[0]/4.0;
    sigma[1] = sigma_[0]/4.0;
    rhot     = rho_[0];
    sigmat   = sigma_[0];    
  }

  if(rhot < MIN_DENS) return;

  /* some handy functions of the total density */
  rhot13   = POW(rhot, 1.0/3.0);
  rhot43   = rhot*rhot13;
  rho83[0] = POW(rho[0], 8.0/3.0);
  rho83[1] = POW(rho[1], 8.0/3.0);

  ZZ     = rhot13/(rhot13 + dd);
  delta  = (cc + dd*ZZ)/rhot13;
  omega  = exp(-cc/rhot13) * ZZ * POW(rhot, -11.0/3.0);

  /* and their derivatives */
  dZZdr    = dd*ZZ*ZZ/(3.0*rhot43);
  ddeltadr = -(cc + dd*ZZ - dd*dd*ZZ*ZZ/rhot13)/(3.0*rhot43);
  domegadr = omega*ZZ*(cc*dd + (cc-10.0*dd)*rhot13 - 11.0*rhot13*rhot13)/(3.0*rhot43*rhot13);

  /* we now add the terms of the xc energy one by one */

  /* t1 */
  {
    FLOAT aux1;

    aux1 = -4.0*AA/rhot;

    t1 = aux1*rho[0]*rho[1]*ZZ;

    if(vrho != NULL){
      for(is=0; is<p->nspin; is++){
	int js = (is==0) ? 1 : 0;

	vrho[is] = aux1*(rho[js]*ZZ - rho[0]*rho[1]*(ZZ/rhot - dZZdr));
      }
    }
  }

  
  /* t2 */
  {
    FLOAT aux1, aux2, aux3, aux4;

    aux1 = -AA*BB;
    aux2 = rho[0]*rho[1];
    aux3 = (47.0 - 7.0*delta)/18.0;
    aux4 = 2.0/3.0*rhot*rhot;

    t2 = aux1*omega*sigmat*(aux2*aux3 - aux4);

    if(vrho != NULL){
      for(is=0; is<p->nspin; is++){
	int js = (is==0) ? 1 : 0;

	vrho[is] += aux1*
	  (domegadr*sigmat*(aux2*aux3 - aux4) +
	   omega*sigmat*(rho[js]*aux3 - aux2*7.0/18.0*ddeltadr - 4.0/3.0*rhot));
      }

      {
	FLOAT aux5 = aux1*omega*(aux2*aux3 - aux4);

	vsigma[0] += aux5/sfact;
	if(p->nspin == XC_POLARIZED){
	  vsigma[1] += 2.0*aux5;
	  vsigma[2] +=     aux5;
	}else{
	  vsigma[0] += aux5/2.0;
	}
      }
    }
  }


  /* t3 */
  {
    FLOAT aux1, aux2;
    aux1 = -AA*BB;
    aux2 = ee*(rho83[0] + rho83[1]);

    t3 = aux1*omega*rho[0]*rho[1]*aux2;

    if(vrho != NULL){
      for(is=0; is<p->nspin; is++){
	int js = (is==0) ? 1 : 0;

	vrho[is] += aux1*domegadr*rho[0]*rho[1]*aux2;
	vrho[is] += aux1* omega * ee
	  * rho[js]*(11.0/3.0*rho83[is] + rho83[js]);
      }
    }
  }


  /* t4 */
  {
    FLOAT aux1, aux2, aux3, aux4;

    aux1 = AA*BB;
    aux2 = rho[0]*rho[1];
    aux3 = 5.0/2.0 - delta/18.0;
    aux4 = sigma[0] + sigma[1];

    t4 = aux1*omega*aux2*aux3*aux4;

    if(vrho != NULL){
      for(is=0; is<p->nspin; is++){
	int js = (is==0) ? 1 : 0;
      
	vrho[is] += aux1*aux4*
	  (domegadr*aux2*aux3 + omega*rho[js]*aux3 -
	   omega*aux2*ddeltadr/18.0);
      }

      {
	FLOAT aux5 = aux1*omega*aux2*aux3;

	vsigma[0] += aux5/sfact;
	if(p->nspin == XC_POLARIZED)
	  vsigma[2] += aux5;
      }
    }
  }


  /* t5 */
  {
    FLOAT aux1, aux2, aux3, aux4;

    aux1 = AA*BB;
    aux2 = rho[0]*rho[1]/(9.0*rhot);
    aux3 = delta - 11.0;
    aux4 = rho[0]*sigma[0] + rho[1]*sigma[1];

    t5 = aux1*omega*aux2*aux3*aux4;

    if(vrho != NULL){
      for(is=0; is<p->nspin; is++){
	int js = (is==0) ? 1 : 0;

	vrho[is] += aux1*
	  ((domegadr*aux2 + omega*(rho[js]*rhot - rho[0]*rho[1])/(9.0*rhot*rhot))*aux3*aux4 +
	   omega*aux2*(ddeltadr*aux4 + aux3*sigma[is]));

	vsigma[is==0 ? 0 : 2] += aux1*omega*aux2*aux3*rho[is]/sfact;
      }
    }
  }

  
  /* t6 */
  {
    FLOAT aux1, aux2, aux3;

    aux1 = -AA*BB;
    aux2 = 2.0/3.0*rhot*rhot*(sigma[0] + sigma[1]);
    aux3 = rho[0]*rho[0]*sigma[1] + rho[1]*rho[1]*sigma[0];

    t6 = aux1*omega*(aux2 - aux3);

    if(vrho != NULL){
      for(is=0; is<p->nspin; is++){
	int js = (is==0) ? 1 : 0;

	vrho[is] += aux1*
	  (domegadr*(aux2 - aux3) + 
	   omega*(4.0/3.0*rhot*(sigma[0] + sigma[1]) - 2*rho[is]*sigma[js]));

	vsigma[is==0 ? 0 : 2] += aux1*omega*(2.0/3.0*rhot*rhot - rho[js]*rho[js])/sfact;
      }
    }
  }

  /* we add all contributions to the total energy */
  *e = (t1 + t2 + t3 + t4 + t5 + t6)/rhot;
}

static void
gga_c_lyp(const void *p_, int np, const FLOAT *rho, const FLOAT *sigma,
	  FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	  FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  int ip;
  const XC(gga_type) *p = p_;

  for(ip=0; ip<np; ip++){
    my_gga_c_lyp(p_, rho, sigma, zk, vrho, vsigma);

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

const XC(func_info_type) XC(func_info_gga_c_lyp) = {
  XC_GGA_C_LYP,
  XC_CORRELATION,
  "Lee, Yang & Parr",
  XC_FAMILY_GGA,
  "C Lee, W Yang and RG Parr, Phys. Rev. B 37, 785 (1988)\n"
  "B Miehlich, A Savin, H Stoll and H Preuss, Chem. Phys. Lett. 157, 200 (1989)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  XC(gga_c_lyp_init), 
  NULL,
  NULL,
  gga_c_lyp
};
