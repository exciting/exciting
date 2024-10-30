/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_gga.c"
#include "funcs_hyb_gga.c"

/* macro to check is a buffer exists */
#define check_out_var(VAR) if(out->VAR == NULL){fprintf(stderr, "error: output variable, out->" #VAR ", is a null pointer\n"); exit(1);}

void
xc_gga_sanity_check(const xc_func_info_type *info, int order, xc_gga_out_params *out)
{
  /* sanity check */
  if(order < 0 || order > 4){
    fprintf(stderr, "Order of derivatives '%d' not implemented\n",
	    order);
    exit(1);
  }
  
  /* sanity check */
  if(out->zk != NULL && !(info->flags & XC_FLAGS_HAVE_EXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of Exc\n",
	    info->name);
    exit(1);
  }

  if(out->vrho != NULL){
    if(!(info->flags & XC_FLAGS_HAVE_VXC)){
      fprintf(stderr, "Functional '%s' does not provide an implementation of vxc\n",
              info->name);
      exit(1);
    }
    check_out_var(vsigma);
  }

  if(out->v2rho2 != NULL){
    if(!(info->flags & XC_FLAGS_HAVE_FXC)){
      fprintf(stderr, "Functional '%s' does not provide an implementation of fxc\n",
              info->name);
      exit(1);
    }
    check_out_var(v2rhosigma); 
    check_out_var(v2sigma2);
  }

  if(out->v3rho3){
    if(!(info->flags & XC_FLAGS_HAVE_KXC)){
      fprintf(stderr, "Functional '%s' does not provide an implementation of kxc\n",
              info->name);
      exit(1);
    }
    check_out_var(v3rho2sigma);
    check_out_var(v3rhosigma2);
    check_out_var(v3sigma3);
  }

  if(out->v4rho4 != NULL){
    if(!(info->flags & XC_FLAGS_HAVE_LXC)){
      fprintf(stderr, "Functional '%s' does not provide an implementation of lxc\n",
              info->name);
      exit(1);
    }
    check_out_var(v4rho3sigma);
    check_out_var(v4rho2sigma2);
    check_out_var(v4rhosigma3);
    check_out_var(v4sigma4);
  }
}

void
xc_gga_initalize(const xc_func_type *func, size_t np, xc_gga_out_params *out)
{
  const xc_dimensions *dim = &(func->dim);

    /* initialize output to zero */
  if(out->zk != NULL)
    libxc_memset(out->zk, 0, dim->zk*np*sizeof(double));

  if(out->vrho != NULL){
    libxc_memset(out->vrho,   0, dim->vrho  *np*sizeof(double));
    libxc_memset(out->vsigma, 0, dim->vsigma*np*sizeof(double));
  }

  if(out->v2rho2 != NULL){
    libxc_memset(out->v2rho2,     0, dim->v2rho2    *np*sizeof(double));
    libxc_memset(out->v2rhosigma, 0, dim->v2rhosigma*np*sizeof(double));
    libxc_memset(out->v2sigma2,   0, dim->v2sigma2  *np*sizeof(double));
  }

  if(out->v3rho3 != NULL){
    libxc_memset(out->v3rho3,      0, dim->v3rho3     *np*sizeof(double));
    libxc_memset(out->v3rho2sigma, 0, dim->v3rho2sigma*np*sizeof(double));
    libxc_memset(out->v3rhosigma2, 0, dim->v3rhosigma2*np*sizeof(double));
    libxc_memset(out->v3sigma3,    0, dim->v3sigma3   *np*sizeof(double));
  }

  if(out->v4rho4 != NULL){
    libxc_memset(out->v4rho4,       0, dim->v4rho4      *np*sizeof(double));
    libxc_memset(out->v4rho3sigma,  0, dim->v4rho3sigma *np*sizeof(double));
    libxc_memset(out->v4rho2sigma2, 0, dim->v4rho2sigma2*np*sizeof(double));
    libxc_memset(out->v4rhosigma3,  0, dim->v4rhosigma3 *np*sizeof(double));
    libxc_memset(out->v4sigma4,     0, dim->v4sigma4    *np*sizeof(double));
   }

}

/* Some useful formulas:

   sigma_st          = grad rho_s . grad rho_t
   zk                = energy density per unit particle

   vrho_s            = d zk / d rho_s
   vsigma_st         = d n*zk / d sigma_st

   v2rho2_st         = d^2 n*zk / d rho_s d rho_t
   v2rhosigma_svx    = d^2 n*zk / d rho_s d sigma_tv
   v2sigma2_stvx     = d^2 n*zk / d sigma_st d sigma_vx

   v3rho3_stv        = d^3 n*zk / d rho_s d rho_t d rho_v
   v3rho2sigma_stvx  = d^3 n*zk / d rho_s d rho_t d sigma_vx
   v3rhosigma2_svxyz = d^3 n*zk / d rho_s d sigma_vx d sigma_yz
   v3sigma3_stvxyz   = d^3 n*zk / d sigma_st d sigma_vx d sigma_yz

if nspin == 2
   rho(2)          = (u, d)
   sigma(3)        = (uu, ud, dd)

   vrho(2)         = (u, d)
   vsigma(3)       = (uu, ud, dd)

   v2rho2(3)       = (u_u, u_d, d_d)
   v2rhosigma(6)   = (u_uu, u_ud, u_dd, d_uu, d_ud, d_dd)
   v2sigma2(6)     = (uu_uu, uu_ud, uu_dd, ud_ud, ud_dd, dd_dd)

   v3rho3(4)       = (u_u_u, u_u_d, u_d_d, d_d_d)
   v3rho2sigma(9)  = (u_u_uu, u_u_ud, u_u_dd, u_d_uu, u_d_ud, u_d_dd, d_d_uu, d_d_ud, d_d_dd)
   v3rhosigma2(12) = (u_uu_uu, u_uu_ud, u_uu_dd, u_ud_ud, u_ud_dd, u_dd_dd, d_uu_uu, d_uu_ud, d_uu_dd, d_ud_ud, d_ud_dd, d_dd_dd)
   v3sigma(10)     = (uu_uu_uu, uu_uu_ud, uu_uu_dd, uu_ud_ud, uu_ud_dd, uu_dd_dd, ud_ud_ud, ud_ud_dd, ud_dd_dd, dd_dd_dd)

*/


void xc_gga_new(const xc_func_type *func, int order, size_t np, const double *rho, const double *sigma,
            xc_gga_out_params *out)
{

  xc_gga_sanity_check(func->info, order, out);
  xc_gga_initalize(func, np, out);
  
  /* call the GGA routines */
  if(func->info->gga != NULL){
    if(func->nspin == XC_UNPOLARIZED){
      if(func->info->gga->unpol[order] != NULL)
        func->info->gga->unpol[order](func, np, rho, sigma, out);
    }else{
      if(func->info->gga->pol[order] != NULL)
        func->info->gga->pol[order](func, np, rho, sigma, out);
    }
  }

  if(func->mix_coef != NULL)
    xc_mix_func(func, np, rho, sigma, NULL, NULL, out->zk, out->vrho, out->vsigma, NULL, NULL,
                out->v2rho2, out->v2rhosigma, NULL, NULL, out->v2sigma2, NULL, NULL, NULL, NULL, NULL,
                out->v3rho3, out->v3rho2sigma, NULL, NULL, out->v3rhosigma2, NULL, NULL, NULL, NULL, NULL,
                out->v3sigma3, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                out->v4rho4, out->v4rho3sigma, NULL, NULL, out->v4rho2sigma2, NULL, NULL, NULL, NULL, NULL,
                out->v4rhosigma3, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                out->v4sigma4, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL);
}

/* old API */
void
xc_gga(const xc_func_type *p, size_t np, const double *rho, const double *sigma,
       double *zk,
       double *vrho, double *vsigma,
       double *v2rho2, double *v2rhosigma, double *v2sigma2,
       double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3,
       double *v4rho4, double *v4rho3sigma, double *v4rho2sigma2, double *v4rhosigma3, double *v4sigma4)
{
  int order = -1;

  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;
  if(v3rho3 != NULL) order = 3;
  if(v4rho4 != NULL) order = 4;

  if(order < 0) return;

  xc_gga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_gga_out_params));
  out.zk     = zk;
  out.vrho   = vrho; out.vsigma = vsigma;
  out.v2rho2 = v2rho2; out.v2rhosigma = v2rhosigma; out.v2sigma2 = v2sigma2;
  out.v3rho3 = v3rho3; out.v3rho2sigma = v3rho2sigma; out.v3rhosigma2 = v3rhosigma2; out.v3sigma3 = v3sigma3;

  out.v4rho4 = v4rho4; out.v4rho3sigma = v4rho3sigma; out.v4rho2sigma2 = v4rho2sigma2; out.v4rhosigma3 = v4rhosigma3; out.v4sigma4 = v4sigma4;

  xc_gga_new(p, order, np, rho, sigma, &out);
}


/* specializations */
void
xc_gga_exc(const xc_func_type *p, size_t np, const double *rho, const double *sigma,
	    double *zk)
{
  xc_gga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_gga_out_params));
  out.zk   = zk;
  
  xc_gga_new(p, 0, np, rho, sigma, &out);
}

void
xc_gga_exc_vxc(const xc_func_type *p, size_t np, const double *rho, const double *sigma,
		double *zk, double *vrho, double *vsigma)
{
  xc_gga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_gga_out_params));
  out.zk   = zk;
  out.vrho = vrho; out.vsigma = vsigma;
  
  xc_gga_new(p, 1, np, rho, sigma, &out);
}

void
xc_gga_exc_vxc_fxc (const xc_func_type *p, size_t np, const double *rho, const double *sigma,
                    double *zk, double *vrho, double *vsigma,
                    double *v2rho2, double *v2rhosigma, double *v2sigma2)
{
  xc_gga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_gga_out_params));
  out.zk     = zk;
  out.vrho   = vrho; out.vsigma = vsigma;
  out.v2rho2 = v2rho2; out.v2rhosigma = v2rhosigma; out.v2sigma2 = v2sigma2;
  
  xc_gga_new(p, 2, np, rho, sigma, &out);
}

void
xc_gga_vxc_fxc (const xc_func_type *p, size_t np, const double *rho, const double *sigma,
                double *vrho, double *vsigma,
                double *v2rho2, double *v2rhosigma, double *v2sigma2)
{
  xc_gga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_gga_out_params));
  out.vrho   = vrho; out.vsigma = vsigma;
  out.v2rho2 = v2rho2; out.v2rhosigma = v2rhosigma; out.v2sigma2 = v2sigma2;
  
  xc_gga_new(p, 2, np, rho, sigma, &out);
}

void
xc_gga_exc_vxc_fxc_kxc (const xc_func_type *p, size_t np, const double *rho, const double *sigma,
                        double *zk, double *vrho, double *vsigma, double *v2rho2, double *v2rhosigma, double *v2sigma2,
                        double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3)
{
  xc_gga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_gga_out_params));
  out.zk     = zk;
  out.vrho   = vrho; out.vsigma = vsigma;
  out.v2rho2 = v2rho2; out.v2rhosigma = v2rhosigma; out.v2sigma2 = v2sigma2;
  out.v3rho3 = v3rho3; out.v3rho2sigma = v3rho2sigma; out.v3rhosigma2 = v3rhosigma2; out.v3sigma3 = v3sigma3;
  
  xc_gga_new(p, 3, np, rho, sigma, &out);
}

void
xc_gga_vxc_fxc_kxc (const xc_func_type *p, size_t np, const double *rho, const double *sigma,
                    double *vrho, double *vsigma, double *v2rho2, double *v2rhosigma, double *v2sigma2,
                    double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3)
{
  xc_gga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_gga_out_params));
  out.vrho   = vrho; out.vsigma = vsigma;
  out.v2rho2 = v2rho2; out.v2rhosigma = v2rhosigma; out.v2sigma2 = v2sigma2;
  out.v3rho3 = v3rho3; out.v3rho2sigma = v3rho2sigma; out.v3rhosigma2 = v3rhosigma2; out.v3sigma3 = v3sigma3;

  xc_gga_new(p, 3, np, rho, sigma, &out);
}

void
xc_gga_vxc(const xc_func_type *p, size_t np, const double *rho, const double *sigma,
           double *vrho, double *vsigma)
{
  xc_gga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_gga_out_params));
  out.vrho   = vrho; out.vsigma = vsigma;

  xc_gga_new(p, 1, np, rho, sigma, &out);
}

void
xc_gga_fxc(const xc_func_type *p, size_t np, const double *rho, const double *sigma,
           double *v2rho2, double *v2rhosigma, double *v2sigma2)
{
  xc_gga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_gga_out_params));
  out.v2rho2 = v2rho2; out.v2rhosigma = v2rhosigma; out.v2sigma2 = v2sigma2;

  xc_gga_new(p, 2, np, rho, sigma, &out);
}

void
xc_gga_kxc(const xc_func_type *p, size_t np, const double *rho, const double *sigma,
           double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3)
{
  xc_gga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_gga_out_params));
  out.v3rho3 = v3rho3; out.v3rho2sigma = v3rho2sigma; out.v3rhosigma2 = v3rhosigma2; out.v3sigma3 = v3sigma3;

  xc_gga_new(p, 3, np, rho, sigma, &out);
}


void
xc_gga_lxc(const xc_func_type *p, size_t np, const double *rho, const double *sigma,
           double *v4rho4, double *v4rho3sigma, double *v4rho2sigma2, double *v4rhosigma3, double *v4sigma4)
{
  xc_gga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_gga_out_params));
  out.v4rho4 = v4rho4; out.v4rho3sigma = v4rho3sigma; out.v4rho2sigma2 = v4rho2sigma2; out.v4rhosigma3 = v4rhosigma3; out.v4sigma4 = v4sigma4;

  xc_gga_new(p, 4, np, rho, sigma, &out);
}
