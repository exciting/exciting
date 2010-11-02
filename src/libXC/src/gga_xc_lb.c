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

/* Note: Do not forget to add a correlation (LDA) functional to the LB94 */
#define XC_GGA_XC_LB 160 /* van Leeuwen & Baerends */

typedef struct{
  int    modified; /* shall we use a modified version */
  FLOAT threshold; /* when to start using the analytic form */
  FLOAT ip;        /* ionization potential of the species */
  FLOAT qtot;      /* total charge in the region */

  FLOAT alpha;     /* the parameters of LB94 */
  FLOAT gamm;
} XC(gga_xc_lb_params);

/************************************************************************
  Calculates van Leeuwen Baerends functional
************************************************************************/

static void
gga_lb_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  assert(p->params == NULL);

  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_LDA_X, p->nspin);

  p->params = malloc(sizeof(XC(gga_xc_lb_params)));
  XC(gga_lb_set_params_)(p, 0, 0.0, 0.0, 0.0);
}


static void
gga_lb_end(void *p_)
{
  XC(gga_type) *p = p_;

  free(p->func_aux[0]);
  free(p->func_aux);
}


void
XC(gga_lb_set_params)(XC(func_type) *p, int modified, FLOAT threshold, FLOAT ip, FLOAT qtot)
{
  assert(p!=NULL && p->gga!=NULL);
  XC(gga_lb_set_params_)(p->gga, modified, threshold, ip, qtot);
}

void
XC(gga_lb_set_params_)(XC(gga_type) *p, int modified, FLOAT threshold, FLOAT ip, FLOAT qtot)
{
  XC(gga_xc_lb_params) *params;

  fflush(stdout);
  assert(p->params != NULL);
  params = (XC(gga_xc_lb_params) *) (p->params);

  params->modified  = modified;
  params->threshold = threshold;
  params->ip        = ip;
  params->qtot      = qtot;

  if(params->modified){
    params->alpha = (params->ip > 0.0) ? 2.0*sqrt(2.0*params->ip) : 0.5;
    params->gamm  = POW(params->qtot, 1.0/3.0)/(2.0*params->alpha);
  }else{
    params->alpha = 0.5;
    params->gamm  = 1.0;
  }
}


void 
XC(gga_lb_modified)(const XC(gga_type) *func, int np, const FLOAT *rho, const FLOAT *sigma, FLOAT r, FLOAT *vrho)
{
  int ip, is;
  FLOAT gdm, x;

  XC(gga_xc_lb_params) *params;

  static const FLOAT beta = 0.05;

  assert(func != NULL);

  assert(func->params != NULL);
  params = (XC(gga_xc_lb_params) *) (func->params);

  XC(lda_vxc)(func->func_aux[0], np, rho, vrho);

  for(ip=0; ip<np; ip++){
    for(is=0; is<func->nspin; is++){
      gdm = sqrt(sigma[(is==0) ? 0 : 2]);

      if(params->modified == 0 || 
	 (rho[is] > params->threshold && gdm > params->threshold)){
	FLOAT f;
      
	if(rho[is] <= MIN_DENS) continue;
	
	x =  gdm/POW(rho[is], 4.0/3.0);
	
	if(x < 300.0) /* the actual functional */	   
	  f = -beta*x*x/(1.0 + 3.0*beta*x*asinh(params->gamm*x));
	else          /* asymptotic expansion */
	  f = -x/(3.0*log(2.0*params->gamm*x));

	vrho[is] += f * POW(rho[is], 1.0/3.0);
	
      }else if(r > 0.0){
	/* the aymptotic expansion of LB94 */
	x = r + (3.0/params->alpha)*
	  log(2.0*params->gamm * params->alpha * POW(params->qtot, -1.0/3.0));
	
	/* x = x + POW(qtot*exp(-alpha*r), 1.0/3.0)/(beta*alpha*alpha); */
	
	vrho[is] -= 1.0/x;
      }
    }
    /* increment pointers */
    rho   += func->n_rho;
    sigma += func->n_sigma;
    
    if(vrho != NULL)
      vrho   += func->n_vrho;

  } /* ip loop */
}


static void 
gga_xc_lb(const void *p_, int np, const FLOAT *rho, const FLOAT *sigma,
	  FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	  FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  XC(gga_lb_modified)((XC(gga_type) *)p_, np, rho, sigma, 0.0, vrho);
}


XC(func_info_type) XC(func_info_gga_xc_lb) = {
  XC_GGA_XC_LB,
  XC_EXCHANGE_CORRELATION,
  "van Leeuwen & Baerends",
  XC_FAMILY_GGA,
  "R van Leeuwen and EJ Baerends, Phys. Rev. A. 49, 2421 (1994)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_VXC,
  gga_lb_init,
  gga_lb_end,
  NULL,
  gga_xc_lb
};

