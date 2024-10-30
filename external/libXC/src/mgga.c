/*
 Copyright (C) 2006-2021 M.A.L. Marques
               2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_mgga.c"
#include "funcs_hyb_mgga.c"

/* macro to check is a buffer exists */
#define check_out_var(VAR) if(out->VAR == NULL){fprintf(stderr, "error: output variable, out->" #VAR ", is a null pointer\n"); exit(1);}

void
xc_mgga_sanity_check(const xc_func_info_type *info, int order, xc_mgga_out_params *out)
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
    if(info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      check_out_var(vlapl);
    }
    if(info->flags & XC_FLAGS_NEEDS_TAU){
      check_out_var(vtau);
    }
  }

  if(out->v2rho2 != NULL){
    if(!(info->flags & XC_FLAGS_HAVE_FXC)){
      fprintf(stderr, "Functional '%s' does not provide an implementation of fxc\n",
              info->name);
      exit(1);
    }
    check_out_var(v2rhosigma); 
    check_out_var(v2sigma2);
    if(info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      check_out_var(v2rholapl);
      check_out_var(v2sigmalapl);
      check_out_var(v2lapl2);
    }
    if(info->flags & XC_FLAGS_NEEDS_TAU){
      check_out_var(v2rhotau);
      check_out_var(v2sigmatau);
      check_out_var(v2tau2);
    }
    if((info->flags & XC_FLAGS_NEEDS_LAPLACIAN) && (info->flags & XC_FLAGS_NEEDS_TAU)){
      check_out_var(v2lapltau);
    }
  }

  if(out->v3rho3 != NULL){
    if(!(info->flags & XC_FLAGS_HAVE_KXC)){
      fprintf(stderr, "Functional '%s' does not provide an implementation of kxc\n",
              info->name);
      exit(1);
    }
    check_out_var(v3rho2sigma);
    check_out_var(v3rhosigma2);
    check_out_var(v3sigma3);
    if(info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      check_out_var(v3rho2lapl);
      check_out_var(v3rhosigmalapl);
      check_out_var(v3rholapl2);
      check_out_var(v3sigma2lapl);
      check_out_var(v3sigmalapl2);
      check_out_var(v3lapl3);
    }
    if(info->flags & XC_FLAGS_NEEDS_TAU){
      check_out_var(v3rho2tau);
      check_out_var(v3rhosigmatau);
      check_out_var(v3rhotau2);
      check_out_var(v3sigma2tau);
      check_out_var(v3sigmatau2);
      check_out_var(v3tau3);
    }
    if((info->flags & XC_FLAGS_NEEDS_LAPLACIAN) && (info->flags & XC_FLAGS_NEEDS_TAU)){
      check_out_var(v3rholapltau);
      check_out_var(v3sigmalapltau);
      check_out_var(v3lapl2tau);
      check_out_var(v3lapltau2);
    }
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
    if(info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      check_out_var(v4rho3lapl);
      check_out_var(v4rho2sigmalapl);
      check_out_var(v4rho2lapl2);
      check_out_var(v4rhosigma2lapl);
      check_out_var(v4rhosigmalapl2);
      check_out_var(v4rholapl3);
      check_out_var(v4sigma3lapl);
      check_out_var(v4sigma2lapl2);
      check_out_var(v4sigmalapl3);
      check_out_var(v4lapl4);
    }
    if(info->flags & XC_FLAGS_NEEDS_TAU){
      check_out_var(v4rho3tau);
      check_out_var(v4rho2sigmatau);
      check_out_var(v4rho2tau2);
      check_out_var(v4rhosigma2tau);
      check_out_var(v4rhosigmatau2);
      check_out_var(v4rhotau3);
      check_out_var(v4sigma3tau);
      check_out_var(v4sigma2tau2);
      check_out_var(v4sigmatau3);
      check_out_var(v4tau4);
    }
    if((info->flags & XC_FLAGS_NEEDS_LAPLACIAN) && (info->flags & XC_FLAGS_NEEDS_TAU)){
      check_out_var(v4rho2lapltau);
      check_out_var(v4rhosigmalapltau);
      check_out_var(v4rholapl2tau);
      check_out_var(v4rholapltau2);
      check_out_var(v4sigma2lapltau);
      check_out_var(v4sigmalapl2tau);
      check_out_var(v4sigmalapltau2);
      check_out_var(v4lapl3tau);
      check_out_var(v4lapl2tau2);
      check_out_var(v4lapltau3); 
    }
  }
}

void
xc_mgga_initalize(const xc_func_type *func, size_t np, xc_mgga_out_params *out)
{
  const xc_dimensions *dim = &(func->dim);

  /* initialize output to zero */
  if(out->zk != NULL)
    libxc_memset(out->zk, 0, dim->zk*np*sizeof(double));

  if(out->vrho != NULL){
    libxc_memset(out->vrho,   0, dim->vrho  *np*sizeof(double));
    libxc_memset(out->vsigma, 0, dim->vsigma*np*sizeof(double));

    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) {
      libxc_memset(out->vlapl,  0, dim->vlapl *np*sizeof(double));
    }
    if(func->info->flags & XC_FLAGS_NEEDS_TAU) {
      libxc_memset(out->vtau,   0, dim->vtau  *np*sizeof(double));
    }
  }

  if(out->v2rho2 != NULL){
    libxc_memset(out->v2rho2,     0, dim->v2rho2     *np*sizeof(double));
    libxc_memset(out->v2rhosigma, 0, dim->v2rhosigma *np*sizeof(double));
    libxc_memset(out->v2sigma2,   0, dim->v2sigma2   *np*sizeof(double));

    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      libxc_memset(out->v2rholapl,   0, dim->v2rholapl  *np*sizeof(double));
      libxc_memset(out->v2sigmalapl, 0, dim->v2sigmalapl*np*sizeof(double));
      libxc_memset(out->v2lapl2,     0, dim->v2lapl2    *np*sizeof(double));
    }
    
    if(func->info->flags & XC_FLAGS_NEEDS_TAU){
      libxc_memset(out->v2rhotau,   0, dim->v2rhotau   *np*sizeof(double));
      libxc_memset(out->v2sigmatau, 0, dim->v2sigmatau *np*sizeof(double));
      libxc_memset(out->v2tau2,     0, dim->v2tau2     *np*sizeof(double));
    }
    
    if((func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) && (func->info->flags & XC_FLAGS_NEEDS_TAU)) {
      libxc_memset(out->v2lapltau,   0, dim->v2lapltau  *np*sizeof(double));
    }
  }

  if(out->v3rho3 != NULL){
    libxc_memset(out->v3rho3,        0, dim->v3rho3       *np*sizeof(double));
    libxc_memset(out->v3rho2sigma,   0, dim->v3rho2sigma  *np*sizeof(double));
    libxc_memset(out->v3rhosigma2,   0, dim->v3rhosigma2  *np*sizeof(double));
    libxc_memset(out->v3sigma3,      0, dim->v3sigma3     *np*sizeof(double));

    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      libxc_memset(out->v3rho2lapl,     0, dim->v3rho2lapl    *np*sizeof(double));
      libxc_memset(out->v3rhosigmalapl, 0, dim->v3rhosigmalapl*np*sizeof(double));
      libxc_memset(out->v3rholapl2,     0, dim->v3rholapl2    *np*sizeof(double));
      libxc_memset(out->v3sigma2lapl,   0, dim->v3sigma2lapl  *np*sizeof(double));
      libxc_memset(out->v3sigmalapl2,   0, dim->v3sigmalapl2  *np*sizeof(double));
      libxc_memset(out->v3lapl3,        0, dim->v3lapl3       *np*sizeof(double));
    }

    if(func->info->flags & XC_FLAGS_NEEDS_TAU){
      libxc_memset(out->v3rho2tau,     0, dim->v3rho2tau    *np*sizeof(double));
      libxc_memset(out->v3rhosigmatau, 0, dim->v3rhosigmatau*np*sizeof(double));
      libxc_memset(out->v3rhotau2,     0, dim->v3rhotau2    *np*sizeof(double));
      libxc_memset(out->v3sigma2tau,   0, dim->v3sigma2tau  *np*sizeof(double));
      libxc_memset(out->v3sigmatau2,   0, dim->v3sigmatau2  *np*sizeof(double));
      libxc_memset(out->v3tau3,        0, dim->v3tau3       *np*sizeof(double));
    }

    if((func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) && (func->info->flags & XC_FLAGS_NEEDS_TAU)) {
      libxc_memset(out->v3rholapltau,   0, dim->v3rholapltau  *np*sizeof(double));
      libxc_memset(out->v3sigmalapltau, 0, dim->v3sigmalapltau*np*sizeof(double));
      libxc_memset(out->v3lapl2tau,     0, dim->v3lapl2tau    *np*sizeof(double));
      libxc_memset(out->v3lapltau2,     0, dim->v3lapltau2    *np*sizeof(double));
    }
  }

  if(out->v4rho4 != NULL){
    libxc_memset(out->v4rho4,         0, dim->v4rho4        *np*sizeof(double));
    libxc_memset(out->v4rho3sigma,    0, dim->v4rho3sigma   *np*sizeof(double));
    libxc_memset(out->v4rho2sigma2,   0, dim->v4rho2sigma2  *np*sizeof(double));
    libxc_memset(out->v4rhosigma3,    0, dim->v4rhosigma3   *np*sizeof(double));
    libxc_memset(out->v4sigma4,       0, dim->v4sigma4      *np*sizeof(double));

    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      libxc_memset(out->v4rho3lapl,        0, dim->v4rho3lapl       *np*sizeof(double));
      libxc_memset(out->v4rho2sigmalapl,   0, dim->v4rho2sigmalapl  *np*sizeof(double));
      libxc_memset(out->v4rho2lapl2,       0, dim->v4rho2lapl2      *np*sizeof(double));
      libxc_memset(out->v4rhosigma2lapl,   0, dim->v4rhosigma2lapl  *np*sizeof(double));
      libxc_memset(out->v4rhosigmalapl2,   0, dim->v4rhosigmalapl2  *np*sizeof(double));
      libxc_memset(out->v4rholapl3,        0, dim->v4rholapl3       *np*sizeof(double));
      libxc_memset(out->v4sigma3lapl,      0, dim->v4sigma3lapl     *np*sizeof(double));
      libxc_memset(out->v4sigma2lapl2,     0, dim->v4sigma2lapl2    *np*sizeof(double));
      libxc_memset(out->v4sigmalapl3,      0, dim->v4sigmalapl3     *np*sizeof(double));
      libxc_memset(out->v4lapl4,           0, dim->v4lapl4          *np*sizeof(double));
    }

    if(func->info->flags & XC_FLAGS_NEEDS_TAU){
      libxc_memset(out->v4rho3tau,      0, dim->v4rho3tau     *np*sizeof(double));
      libxc_memset(out->v4rho2sigmatau, 0, dim->v4rho2sigmatau*np*sizeof(double));
      libxc_memset(out->v4rho2tau2,     0, dim->v4rho2tau2    *np*sizeof(double));
      libxc_memset(out->v4rhosigma2tau, 0, dim->v4rhosigma2tau*np*sizeof(double));
      libxc_memset(out->v4rhosigmatau2, 0, dim->v4rhosigmatau2*np*sizeof(double));
      libxc_memset(out->v4rhotau3,      0, dim->v4rhotau3     *np*sizeof(double));
      libxc_memset(out->v4sigma3tau,    0, dim->v4sigma3tau   *np*sizeof(double));
      libxc_memset(out->v4sigma2tau2,   0, dim->v4sigma2tau2  *np*sizeof(double));
      libxc_memset(out->v4sigmatau3,    0, dim->v4sigmatau3   *np*sizeof(double));
      libxc_memset(out->v4tau4,         0, dim->v4tau4        *np*sizeof(double));
    }

    if((func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) && (func->info->flags & XC_FLAGS_NEEDS_TAU)) {
      libxc_memset(out->v4rho2lapltau,     0, dim->v4rho2lapltau    *np*sizeof(double));
      libxc_memset(out->v4rhosigmalapltau, 0, dim->v4rhosigmalapltau*np*sizeof(double));
      libxc_memset(out->v4rholapl2tau,     0, dim->v4rholapl2tau    *np*sizeof(double));
      libxc_memset(out->v4rholapltau2,     0, dim->v4rholapltau2    *np*sizeof(double));
      libxc_memset(out->v4sigma2lapltau,   0, dim->v4sigma2lapltau  *np*sizeof(double));
      libxc_memset(out->v4sigmalapl2tau,   0, dim->v4sigmalapl2tau  *np*sizeof(double));
      libxc_memset(out->v4sigmalapltau2,   0, dim->v4sigmalapltau2  *np*sizeof(double));
      libxc_memset(out->v4lapl3tau,        0, dim->v4lapl3tau       *np*sizeof(double));
      libxc_memset(out->v4lapl2tau2,       0, dim->v4lapl2tau2      *np*sizeof(double));
      libxc_memset(out->v4lapltau3,        0, dim->v4lapltau3       *np*sizeof(double));
    }
  }

}

void xc_mgga_new(const xc_func_type *func, int order, size_t np,
                const double *rho, const double *sigma, const double *lapl, const double *tau,
                xc_mgga_out_params *out)
{
  xc_mgga_sanity_check(func->info, order, out);
  xc_mgga_initalize(func, np, out);

  /* call the mGGA routines */
  if(func->info->mgga != NULL){
    if(func->nspin == XC_UNPOLARIZED){
      if(func->info->mgga->unpol[order] != NULL)
        func->info->mgga->unpol[order](func, np, rho, sigma, lapl, tau, out);
    }else{
      if(func->info->mgga->pol[order] != NULL)
        func->info->mgga->pol[order](func, np, rho, sigma, lapl, tau, out);
    }
  }
    
  if(func->mix_coef != NULL)
    xc_mix_func(func, np, rho, sigma, lapl, tau,
                out->zk, out->vrho, out->vsigma, out->vlapl, out->vtau,
                out->v2rho2, out->v2rhosigma, out->v2rholapl, out->v2rhotau,
                out->v2sigma2, out->v2sigmalapl, out->v2sigmatau, out->v2lapl2,
                out->v2lapltau, out->v2tau2,
                out->v3rho3, out->v3rho2sigma, out->v3rho2lapl, out->v3rho2tau,
                out->v3rhosigma2, out->v3rhosigmalapl, out->v3rhosigmatau,
                out->v3rholapl2, out->v3rholapltau, out->v3rhotau2, out->v3sigma3,
                out->v3sigma2lapl, out->v3sigma2tau, out->v3sigmalapl2, out->v3sigmalapltau,
                out->v3sigmatau2, out->v3lapl3, out->v3lapl2tau, out->v3lapltau2,
                out->v3tau3,
                out->v4rho4, out->v4rho3sigma, out->v4rho3lapl, out->v4rho3tau,
                out->v4rho2sigma2, out->v4rho2sigmalapl, out->v4rho2sigmatau,
                out->v4rho2lapl2, out->v4rho2lapltau, out->v4rho2tau2, out->v4rhosigma3,
                out->v4rhosigma2lapl, out->v4rhosigma2tau, out->v4rhosigmalapl2,
                out->v4rhosigmalapltau, out->v4rhosigmatau2, out->v4rholapl3,
                out->v4rholapl2tau, out->v4rholapltau2, out->v4rhotau3, out->v4sigma4,
                out->v4sigma3lapl, out->v4sigma3tau, out->v4sigma2lapl2, out->v4sigma2lapltau,
                out->v4sigma2tau2, out->v4sigmalapl3, out->v4sigmalapl2tau,
                out->v4sigmalapltau2, out->v4sigmatau3, out->v4lapl4, out->v4lapl3tau,
                out->v4lapl2tau2, out->v4lapltau3, out->v4tau4
                );
}

/* old API */
#define SET_ORDER_0                          \
  out.zk = zk;

#define SET_ORDER_1                          \
  out.vrho = vrho;                           \
  out.vsigma = vsigma;                       \
  out.vlapl = vlapl;                         \
  out.vtau = vtau;

#define SET_ORDER_2                          \
  out.v2rho2 = v2rho2;                       \
  out.v2rhosigma = v2rhosigma;               \
  out.v2rholapl = v2rholapl;                 \
  out.v2rhotau = v2rhotau;                   \
  out.v2sigma2 = v2sigma2;                   \
  out.v2sigmalapl = v2sigmalapl;             \
  out.v2sigmatau = v2sigmatau;               \
  out.v2lapl2 = v2lapl2;                     \
  out.v2lapltau = v2lapltau;                 \
  out.v2tau2 = v2tau2;

#define SET_ORDER_3                          \
  out.v3rho3 = v3rho3;                       \
  out.v3rho2sigma = v3rho2sigma;             \
  out.v3rho2lapl = v3rho2lapl;               \
  out.v3rho2tau = v3rho2tau;                 \
  out.v3rhosigma2 = v3rhosigma2;             \
  out.v3rhosigmalapl = v3rhosigmalapl;       \
  out.v3rhosigmatau = v3rhosigmatau;         \
  out.v3rholapl2 = v3rholapl2;               \
  out.v3rholapltau = v3rholapltau;           \
  out.v3rhotau2 = v3rhotau2;                 \
  out.v3sigma3 = v3sigma3;                   \
  out.v3sigma2lapl = v3sigma2lapl;           \
  out.v3sigma2tau = v3sigma2tau;             \
  out.v3sigmalapl2 = v3sigmalapl2;           \
  out.v3sigmalapltau = v3sigmalapltau;       \
  out.v3sigmatau2 = v3sigmatau2;             \
  out.v3lapl3 = v3lapl3;                     \
  out.v3lapl2tau = v3lapl2tau;               \
  out.v3lapltau2 = v3lapltau2;               \
  out.v3tau3 = v3tau3;

#define SET_ORDER_4                          \
  out.v4rho4 = v4rho4;                       \
  out.v4rho3sigma = v4rho3sigma;             \
  out.v4rho3lapl = v4rho3lapl;               \
  out.v4rho3tau = v4rho3tau;                 \
  out.v4rho2sigma2 = v4rho2sigma2;           \
  out.v4rho2sigmalapl = v4rho2sigmalapl;     \
  out.v4rho2sigmatau = v4rho2sigmatau;       \
  out.v4rho2lapl2 = v4rho2lapl2;             \
  out.v4rho2lapltau = v4rho2lapltau;         \
  out.v4rho2tau2 = v4rho2tau2;               \
  out.v4rhosigma3 = v4rhosigma3;             \
  out.v4rhosigma2lapl = v4rhosigma2lapl;     \
  out.v4rhosigma2tau = v4rhosigma2tau;       \
  out.v4rhosigmalapl2 = v4rhosigmalapl2;     \
  out.v4rhosigmalapltau = v4rhosigmalapltau; \
  out.v4rhosigmatau2 = v4rhosigmatau2;       \
  out.v4rholapl3 = v4rholapl3;               \
  out.v4rholapl2tau = v4rholapl2tau;         \
  out.v4rholapltau2 = v4rholapltau2;         \
  out.v4rhotau3 = v4rhotau3;                 \
  out.v4sigma4 = v4sigma4;                   \
  out.v4sigma3lapl = v4sigma3lapl;           \
  out.v4sigma3tau = v4sigma3tau;             \
  out.v4sigma2lapl2 = v4sigma2lapl2;         \
  out.v4sigma2lapltau = v4sigma2lapltau;     \
  out.v4sigma2tau2 = v4sigma2tau2;           \
  out.v4sigmalapl3 = v4sigmalapl3;           \
  out.v4sigmalapl2tau = v4sigmalapl2tau;     \
  out.v4sigmalapltau2 = v4sigmalapltau2;     \
  out.v4sigmatau3 = v4sigmatau3;             \
  out.v4lapl4 = v4lapl4;                     \
  out.v4lapl3tau = v4lapl3tau;               \
  out.v4lapl2tau2 = v4lapl2tau2;             \
  out.v4lapltau3 = v4lapltau3;               \
  out.v4tau4 = v4tau4;

void
xc_mgga(const xc_func_type *p, size_t np,
        const double *rho, const double *sigma, const double *lapl, const double *tau,
        double *zk,
        double *vrho, double *vsigma, double *vlapl, double *vtau,
        double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau, double *v2sigma2,
        double *v2sigmalapl, double *v2sigmatau, double *v2lapl2, double *v2lapltau, double *v2tau2,
        double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau, double *v3rhosigma2,
        double *v3rhosigmalapl, double *v3rhosigmatau, double *v3rholapl2, double *v3rholapltau,
        double *v3rhotau2, double *v3sigma3, double *v3sigma2lapl, double *v3sigma2tau,
        double *v3sigmalapl2, double *v3sigmalapltau, double *v3sigmatau2, double *v3lapl3,
        double *v3lapl2tau, double *v3lapltau2, double *v3tau3,
        double *v4rho4, double *v4rho3sigma, double *v4rho3lapl, double *v4rho3tau, double *v4rho2sigma2,
        double *v4rho2sigmalapl, double *v4rho2sigmatau, double *v4rho2lapl2, double *v4rho2lapltau,
        double *v4rho2tau2, double *v4rhosigma3, double *v4rhosigma2lapl, double *v4rhosigma2tau,
        double *v4rhosigmalapl2, double *v4rhosigmalapltau, double *v4rhosigmatau2,
        double *v4rholapl3, double *v4rholapl2tau, double *v4rholapltau2, double *v4rhotau3,
        double *v4sigma4, double *v4sigma3lapl, double *v4sigma3tau, double *v4sigma2lapl2,
        double *v4sigma2lapltau, double *v4sigma2tau2, double *v4sigmalapl3, double *v4sigmalapl2tau,
        double *v4sigmalapltau2, double *v4sigmatau3, double *v4lapl4, double *v4lapl3tau,
        double *v4lapl2tau2, double *v4lapltau3, double *v4tau4
        )
{
  int order = -1;

  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;
  if(v3rho3 != NULL) order = 3;
  if(v4rho4 != NULL) order = 4;

  if(order < 0) return;

  xc_mgga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_mgga_out_params));

  SET_ORDER_0;
  SET_ORDER_1;
  SET_ORDER_2;
  SET_ORDER_3;
  SET_ORDER_4;

  /* call new API */
  xc_mgga_new(p, order, np, rho, sigma, lapl, tau, &out);
}


/* specializations */
void
xc_mgga_exc(const xc_func_type *p, size_t np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *zk)
{
  xc_mgga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_mgga_out_params));
  SET_ORDER_0;

  xc_mgga_new(p, 0, np, rho, sigma, lapl, tau, &out);
}

void
xc_mgga_exc_vxc(const xc_func_type *p, size_t np,
                const double *rho, const double *sigma, const double *lapl, const double *tau,
                double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau)
{
  xc_mgga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_mgga_out_params));
  SET_ORDER_0;
  SET_ORDER_1;

  xc_mgga_new(p, 1, np, rho, sigma, lapl, tau, &out);
}

void xc_mgga_exc_vxc_fxc(const xc_func_type *p, size_t np,
                         const double *rho, const double *sigma, const double *lapl, const double *tau,
                         double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau,
                         double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                         double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                         double *v2lapltau, double *v2tau2)
{
  xc_mgga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_mgga_out_params));
  SET_ORDER_0;
  SET_ORDER_1;
  SET_ORDER_2;

  xc_mgga_new(p, 2, np, rho, sigma, lapl, tau, &out);
}

void xc_mgga_vxc_fxc(const xc_func_type *p, size_t np,
                         const double *rho, const double *sigma, const double *lapl, const double *tau,
                         double *vrho, double *vsigma, double *vlapl, double *vtau,
                         double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                         double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                         double *v2lapltau, double *v2tau2)
{
  xc_mgga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_mgga_out_params));
  SET_ORDER_1;
  SET_ORDER_2;

  xc_mgga_new(p, 2, np, rho, sigma, lapl, tau, &out);
}

void xc_mgga_exc_vxc_fxc_kxc(const xc_func_type *p, size_t np,
                             const double *rho, const double *sigma, const double *lapl, const double *tau,
                             double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau,
                             double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                             double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                             double *v2lapltau, double *v2tau2,
                             double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
                             double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
                             double *v3rholapl2, double *v3rholapltau, double *v3rhotau2, double *v3sigma3,
                             double *v3sigma2lapl, double *v3sigma2tau, double *v3sigmalapl2, double *v3sigmalapltau,
                             double *v3sigmatau2, double *v3lapl3, double *v3lapl2tau, double *v3lapltau2,
                             double *v3tau3)
{
  xc_mgga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_mgga_out_params));
  SET_ORDER_0;
  SET_ORDER_1;
  SET_ORDER_2;
  SET_ORDER_3;

  xc_mgga_new(p, 3, np, rho, sigma, lapl, tau, &out);
}

void xc_mgga_vxc_fxc_kxc(const xc_func_type *p, size_t np,
                         const double *rho, const double *sigma, const double *lapl, const double *tau,
                         double *vrho, double *vsigma, double *vlapl, double *vtau,
                         double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                         double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                         double *v2lapltau, double *v2tau2,
                         double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
                         double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
                         double *v3rholapl2, double *v3rholapltau, double *v3rhotau2, double *v3sigma3,
                         double *v3sigma2lapl, double *v3sigma2tau, double *v3sigmalapl2, double *v3sigmalapltau,
                         double *v3sigmatau2, double *v3lapl3, double *v3lapl2tau, double *v3lapltau2,
                         double *v3tau3)
{
  xc_mgga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_mgga_out_params));
  SET_ORDER_1;
  SET_ORDER_2;
  SET_ORDER_3;

  xc_mgga_new(p, 3, np, rho, sigma, lapl, tau, &out);
}


void
xc_mgga_vxc(const xc_func_type *p, size_t np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *vrho, double *vsigma, double *vlapl, double *vtau)
{
  xc_mgga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_mgga_out_params));
  SET_ORDER_1;

  xc_mgga_new(p, 1, np, rho, sigma, lapl, tau, &out);
}

void
xc_mgga_fxc(const xc_func_type *p, size_t np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
            double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
            double *v2lapltau, double *v2tau2)
{
  xc_mgga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_mgga_out_params));
  SET_ORDER_2;

  xc_mgga_new(p, 2, np, rho, sigma, lapl, tau, &out);
}

void xc_mgga_kxc(const xc_func_type *p, size_t np,
                 const double *rho, const double *sigma, const double *lapl, const double *tau,
                 double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
                 double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
                 double *v3rholapl2, double *v3rholapltau,  double *v3rhotau2,
                 double *v3sigma3, double *v3sigma2lapl, double *v3sigma2tau,
                 double *v3sigmalapl2, double *v3sigmalapltau, double *v3sigmatau2,
                 double *v3lapl3, double *v3lapl2tau, double *v3lapltau2, double *v3tau3)
{
  xc_mgga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_mgga_out_params));
  SET_ORDER_3;

  xc_mgga_new(p, 3, np, rho, sigma, lapl, tau, &out);
}

void xc_mgga_lxc(const xc_func_type *p, size_t np,
                 const double *rho, const double *sigma, const double *lapl, const double *tau,
                 double *v4rho4, double *v4rho3sigma, double *v4rho3lapl, double *v4rho3tau, double *v4rho2sigma2,
                 double *v4rho2sigmalapl, double *v4rho2sigmatau, double *v4rho2lapl2, double *v4rho2lapltau,
                 double *v4rho2tau2, double *v4rhosigma3, double *v4rhosigma2lapl, double *v4rhosigma2tau,
                 double *v4rhosigmalapl2, double *v4rhosigmalapltau, double *v4rhosigmatau2,
                 double *v4rholapl3, double *v4rholapl2tau, double *v4rholapltau2, double *v4rhotau3,
                 double *v4sigma4, double *v4sigma3lapl, double *v4sigma3tau, double *v4sigma2lapl2,
                 double *v4sigma2lapltau, double *v4sigma2tau2, double *v4sigmalapl3, double *v4sigmalapl2tau,
                 double *v4sigmalapltau2, double *v4sigmatau3, double *v4lapl4, double *v4lapl3tau,
                 double *v4lapl2tau2, double *v4lapltau3, double *v4tau4)
{
  xc_mgga_out_params out;
  libxc_memset(&out, 0, sizeof(xc_mgga_out_params));
  SET_ORDER_4;

  xc_mgga_new(p, 4, np, rho, sigma, lapl, tau, &out);
}
