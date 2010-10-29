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
#include <string.h>
#include <assert.h>

#include "util.h"
#include "funcs_mgga.c"

/* initialization */
int XC(mgga_init)(XC(func_type) *p, const XC(func_info_type) *info, int nspin)
{
  XC(mgga_type) *func;

  assert(p != NULL && p->mgga != NULL);
  func = p->mgga;

  /* initialize structure */
  func->info   = info;
  func->nspin  = nspin;
  func->params = NULL;
  func->func   = 0;

  func->n_func_aux = 0;
  func->func_aux   = NULL;
  func->mix_coef   = NULL;

  /* initialize spin counters */
  func->n_zk  = 1;
  func->n_rho = func->n_vrho = func->nspin;
  func->n_tau = func->n_vtau = func->nspin;
  func->n_lapl_rho = func->n_vlapl_rho = func->nspin;
  if(func->nspin == XC_UNPOLARIZED){
    func->n_sigma  = func->n_vsigma = 1;
    func->n_v2rho2 = func->n_v2rhosigma = func->n_v2sigma2 = 1;
  }else{
    func->n_sigma      = func->n_vsigma = func->n_v2rho2 = 3;
    func->n_v2rhosigma = func->n_v2sigma2 = 6;
  }

  /* see if we need to initialize the functional */
  if(func->info->init != NULL)
    func->info->init(func);
  return 0;
}


void XC(mgga_end)(XC(func_type) *p)
{
  XC(mgga_type) *func;

  assert(p != NULL && p->mgga != NULL);
  func = p->mgga;

  /* call internal termination routine */
  if(func->info->end != NULL)
    func->info->end(func);

  /* terminate any auxiliary functional */
  if(func->n_func_aux > 0){
    int ii;
    
    for(ii=0; ii<func->n_func_aux; ii++){
      XC(func_end)(func->func_aux[ii]);
      free(func->func_aux[ii]);
    }
    free(func->func_aux);
  }

  if(func->mix_coef != NULL){
    free(func->mix_coef);
    func->mix_coef = NULL;
  }

  /* deallocate any used parameter */
  if(func->params != NULL){
    free(func->params);
    func->params = NULL;
  }
}


void 
XC(mgga)(const XC(func_type) *p, int np,
	 const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
	 FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau,
	 FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2)
{
  XC(mgga_type) *func;

  assert(p != NULL && p->mgga != NULL);
  func = p->mgga;

  /* sanity check */
  if(zk != NULL && !(func->info->flags & XC_FLAGS_HAVE_EXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of Exc",
	    func->info->name);
    exit(1);
  }

  if(vrho != NULL && !(func->info->flags & XC_FLAGS_HAVE_VXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of vxc",
	    func->info->name);
    exit(1);
  }

  if(v2rho2 != NULL && !(func->info->flags & XC_FLAGS_HAVE_FXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of fxc",
	    func->info->name);
    exit(1);
  }

  /* initialize output to zero */
  if(zk != NULL)
    memset(zk, 0, func->n_zk*np*sizeof(FLOAT));

  if(vrho != NULL){
    assert(vsigma != NULL);

    memset(vrho,      0, func->n_vrho     *np*sizeof(FLOAT));
    memset(vsigma,    0, func->n_vsigma   *np*sizeof(FLOAT));
    memset(vtau,      0, func->n_vtau     *np*sizeof(FLOAT));
    memset(vlapl_rho, 0, func->n_vlapl_rho*np*sizeof(FLOAT));
  }

  if(v2rho2 != NULL){
    /* warning : lapl_rho terms missing here */
    assert(v2rhosigma!=NULL && v2sigma2!=NULL && v2rhotau!=NULL && v2tausigma!=NULL && v2tau2!=NULL);

    memset(v2rho2,     0, func->n_v2rho2    *np*sizeof(FLOAT));
    memset(v2rhosigma, 0, func->n_v2rhosigma*np*sizeof(FLOAT));
    memset(v2sigma2,   0, func->n_v2sigma2  *np*sizeof(FLOAT));
    memset(v2rhotau,   0, func->n_v2rhotau  *np*sizeof(FLOAT));
    memset(v2tausigma, 0, func->n_v2tausigma*np*sizeof(FLOAT));
    memset(v2tau2,     0, func->n_v2tau2    *np*sizeof(FLOAT));
  }

  /* call functional */
  if(func->info->mgga != NULL)
    func->info->mgga(func, np, rho, sigma, lapl_rho, tau, zk, vrho, vsigma, vlapl_rho, vtau, 
		     v2rho2, v2rhosigma, v2sigma2, v2rhotau, v2tausigma, v2tau2);

  /* Mixing still not implemented for mggas
  if(func->mix_coef != NULL){
    XC(mix_func)(p, func->n_func_aux, func->func_aux, func->mix_coef, 
		 np, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);
  }
  */
}

/* especializations */
inline void 
XC(mgga_exc)(const XC(func_type) *p, int np, 
	     const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
	     FLOAT *zk)
{
  XC(mgga)(p, np, rho, sigma, tau, zk, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

inline void 
XC(mgga_exc_vxc)(const XC(func_type) *p, int np,
		 const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
		 FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau)
{
  XC(mgga)(p, np, rho, sigma, lapl_rho, tau, zk, vrho, vsigma, vlapl_rho, vtau, NULL, NULL, NULL, NULL, NULL, NULL);
}

inline void 
XC(mgga_vxc)(const XC(func_type) *p, int np,
	     const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
	     FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau)
{
  XC(mgga)(p, np, rho, sigma, lapl_rho, tau, NULL, vrho, vsigma, vlapl_rho, vtau, NULL, NULL, NULL, NULL, NULL, NULL);
}

inline void 
XC(mgga_fxc)(const XC(func_type) *p, int np,
	     const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
	     FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2)
{
  XC(mgga)(p, np, rho, sigma, lapl_rho, tau, NULL, NULL, NULL, NULL, NULL, v2rho2, v2rhosigma, v2sigma2, v2rhotau, v2tausigma, v2tau2);
}


