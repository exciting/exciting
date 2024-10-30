/*
  Copyright (C) 2006-2007 M.A.L. Marques
                2018-2019 Susi Lehtola
                2019 X. Andrade

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define is_mgga(id)   ((id) == XC_FAMILY_MGGA || (id) == XC_FAMILY_HYB_MGGA)
#define is_gga(id)    ((id) == XC_FAMILY_GGA || (id) == XC_FAMILY_HYB_GGA || is_mgga(id))
#define is_lda(id)    ((id) == XC_FAMILY_LDA || (id) == XC_FAMILY_HYB_LDA ||  is_gga(id))
#define safe_free(pt) if(pt != NULL) libxc_free(pt)

/* Due to the fact that some functionals do not depend on tau or lapl,
   some variables may remain uninitialized when evaluating the ked or
   mgga functionals. Therefore, *all* variables should be explicitly
   initialized to zero here!
*/
void
xc_mgga_vars_allocate_all(int family, size_t np, const xc_dimensions *dim,
                     int do_zk, int do_vrho, int do_v2rho2, int do_v3rho3, int do_v4rho4,
                     double **zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double **, ))
{
  /* allocate buffers */
  if(do_zk){
    *zk = (double *) libxc_malloc(sizeof(double)*np*dim->zk);
    libxc_memset(*zk, 0, dim->zk*np*sizeof(double));
  }

#ifndef XC_DONT_COMPILE_VXC
  if(do_vrho){
    *vrho = (double *) libxc_malloc(sizeof(double)*np*dim->vrho);
    libxc_memset(*vrho,   0, dim->vrho  *np*sizeof(double));
    if(is_gga(family)){
      *vsigma = (double *) libxc_malloc(sizeof(double)*np*dim->vsigma);
      libxc_memset(*vsigma, 0, dim->vsigma*np*sizeof(double));
    }
    if(is_mgga(family)){
      *vlapl = (double *) libxc_malloc(sizeof(double)*np*dim->vlapl);
      libxc_memset(*vlapl,  0, dim->vlapl *np*sizeof(double));
      *vtau  = (double *) libxc_malloc(sizeof(double)*np*dim->vtau);
      libxc_memset(*vtau,   0, dim->vtau  *np*sizeof(double));
    }
  }

#ifndef XC_DONT_COMPILE_FXC
  if(do_v2rho2){
    *v2rho2 = (double *) libxc_malloc(sizeof(double)*np*dim->v2rho2);
     libxc_memset(*v2rho2,     0, dim->v2rho2     *np*sizeof(double));
    if(is_gga(family)){
      *v2rhosigma  = (double *) libxc_malloc(sizeof(double)*np*dim->v2rhosigma);
      libxc_memset(*v2rhosigma, 0, dim->v2rhosigma *np*sizeof(double));
      *v2sigma2    = (double *) libxc_malloc(sizeof(double)*np*dim->v2sigma2);
      libxc_memset(*v2sigma2,   0, dim->v2sigma2   *np*sizeof(double));
    }
    if(is_mgga(family)){
      *v2rholapl   = (double *) libxc_malloc(sizeof(double)*np*dim->v2rholapl);
      libxc_memset(*v2rholapl,   0, dim->v2rholapl  *np*sizeof(double));
      *v2rhotau    = (double *) libxc_malloc(sizeof(double)*np*dim->v2rhotau);
      libxc_memset(*v2rhotau,   0, dim->v2rhotau   *np*sizeof(double));
      *v2sigmalapl = (double *) libxc_malloc(sizeof(double)*np*dim->v2sigmalapl);
      libxc_memset(*v2sigmalapl, 0, dim->v2sigmalapl*np*sizeof(double));
      *v2sigmatau  = (double *) libxc_malloc(sizeof(double)*np*dim->v2sigmatau);
      libxc_memset(*v2sigmatau, 0, dim->v2sigmatau *np*sizeof(double));
      *v2lapl2     = (double *) libxc_malloc(sizeof(double)*np*dim->v2lapl2);
      libxc_memset(*v2lapl2,     0, dim->v2lapl2    *np*sizeof(double));
      *v2lapltau   = (double *) libxc_malloc(sizeof(double)*np*dim->v2lapltau);
      libxc_memset(*v2lapltau,   0, dim->v2lapltau  *np*sizeof(double));
      *v2tau2      = (double *) libxc_malloc(sizeof(double)*np*dim->v2tau2);
      libxc_memset(*v2tau2,     0, dim->v2tau2     *np*sizeof(double));
    }
  }

#ifndef XC_DONT_COMPILE_KXC
  if(do_v3rho3){
    *v3rho3      = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho3);
    libxc_memset(*v3rho3,        0, dim->v3rho3       *np*sizeof(double));
    if(is_gga(family)){
      *v3rho2sigma = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho2sigma);
      libxc_memset(*v3rho2sigma,   0, dim->v3rho2sigma  *np*sizeof(double));
      *v3rhosigma2 = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhosigma2);
      libxc_memset(*v3rhosigma2,   0, dim->v3rhosigma2  *np*sizeof(double));
      *v3sigma3    = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigma3);
      libxc_memset(*v3sigma3,      0, dim->v3sigma3     *np*sizeof(double));
    }
    if(is_mgga(family)){
      *v3rho2lapl     = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho2lapl);
      libxc_memset(*v3rho2lapl,     0, dim->v3rho2lapl    *np*sizeof(double));
      *v3rho2tau      = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho2tau);
      libxc_memset(*v3rho2tau,     0, dim->v3rho2tau    *np*sizeof(double));
      *v3rhosigmalapl = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhosigmalapl);
      libxc_memset(*v3rhosigmalapl, 0, dim->v3rhosigmalapl*np*sizeof(double));
      *v3rhosigmatau  = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhosigmatau);
      libxc_memset(*v3rhosigmatau, 0, dim->v3rhosigmatau*np*sizeof(double));
      *v3rholapl2     = (double *) libxc_malloc(sizeof(double)*np*dim->v3rholapl2);
      libxc_memset(*v3rholapl2,     0, dim->v3rholapl2    *np*sizeof(double));
      *v3rholapltau   = (double *) libxc_malloc(sizeof(double)*np*dim->v3rholapltau);
      libxc_memset(*v3rholapltau,   0, dim->v3rholapltau  *np*sizeof(double));
      *v3rhotau2      = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhotau2);
      libxc_memset(*v3rhotau2,     0, dim->v3rhotau2    *np*sizeof(double));
      *v3sigma2lapl   = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigma2lapl);
      libxc_memset(*v3sigma2lapl,   0, dim->v3sigma2lapl  *np*sizeof(double));
      *v3sigma2tau    = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigma2tau);
      libxc_memset(*v3sigma2tau,   0, dim->v3sigma2tau  *np*sizeof(double));
      *v3sigmalapl2   = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigmalapl2);
      libxc_memset(*v3sigmalapl2,   0, dim->v3sigmalapl2  *np*sizeof(double));
      *v3sigmalapltau = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigmalapltau);
      libxc_memset(*v3sigmalapltau, 0, dim->v3sigmalapltau*np*sizeof(double));
      *v3sigmatau2    = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigmatau2);
      libxc_memset(*v3sigmatau2,   0, dim->v3sigmatau2  *np*sizeof(double));
      *v3lapl3        = (double *) libxc_malloc(sizeof(double)*np*dim->v3lapl3);
      libxc_memset(*v3lapl3,        0, dim->v3lapl3       *np*sizeof(double));
      *v3lapl2tau     = (double *) libxc_malloc(sizeof(double)*np*dim->v3lapl2tau);
      libxc_memset(*v3lapl2tau,     0, dim->v3lapl2tau    *np*sizeof(double));
      *v3lapltau2     = (double *) libxc_malloc(sizeof(double)*np*dim->v3lapltau2);
      libxc_memset(*v3lapltau2,     0, dim->v3lapltau2    *np*sizeof(double));
      *v3tau3         = (double *) libxc_malloc(sizeof(double)*np*dim->v3tau3);
      libxc_memset(*v3tau3,        0, dim->v3tau3       *np*sizeof(double));
    }
  }

#ifndef XC_DONT_COMPILE_LXC
  if(do_v4rho4){
    *v4rho4            = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho4);
    libxc_memset(*v4rho4,         0, dim->v4rho4        *np*sizeof(double));
    if(is_gga(family)){
      *v4rho3sigma       = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho3sigma);
      libxc_memset(*v4rho3sigma,    0, dim->v4rho3sigma   *np*sizeof(double));
      *v4rho2sigma2      = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2sigma2);
      libxc_memset(*v4rho2sigma2,   0, dim->v4rho2sigma2  *np*sizeof(double));
      *v4rhosigma3       = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigma3);
      libxc_memset(*v4rhosigma3,    0, dim->v4rhosigma3   *np*sizeof(double));
      *v4sigma4          = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma4);
      libxc_memset(*v4sigma4,       0, dim->v4sigma4      *np*sizeof(double));
    }
    if(is_mgga(family)){
      *v4rho3lapl        = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho3lapl);
      libxc_memset(*v4rho3lapl,        0, dim->v4rho3lapl       *np*sizeof(double));
      *v4rho3tau         = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho3tau);
      libxc_memset(*v4rho3tau,      0, dim->v4rho3tau     *np*sizeof(double));
      *v4rho2sigmalapl   = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2sigmalapl);
      libxc_memset(*v4rho2sigmalapl,   0, dim->v4rho2sigmalapl  *np*sizeof(double));
      *v4rho2sigmatau    = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2sigmatau);
      libxc_memset(*v4rho2sigmatau, 0, dim->v4rho2sigmatau*np*sizeof(double));
      *v4rho2lapl2       = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2lapl2);
      libxc_memset(*v4rho2lapl2,       0, dim->v4rho2lapl2      *np*sizeof(double));
      *v4rho2lapltau     = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2lapltau);
      libxc_memset(*v4rho2lapltau,     0, dim->v4rho2lapltau    *np*sizeof(double));
      *v4rho2tau2        = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2tau2);
      libxc_memset(*v4rho2tau2,     0, dim->v4rho2tau2    *np*sizeof(double));
      *v4rhosigma2lapl   = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigma2lapl);
      libxc_memset(*v4rhosigma2lapl,   0, dim->v4rhosigma2lapl  *np*sizeof(double));
      *v4rhosigma2tau    = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigma2tau);
      libxc_memset(*v4rho2sigmatau, 0, dim->v4rho2sigmatau*np*sizeof(double));
      *v4rhosigmalapl2   = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigmalapl2);
      libxc_memset(*v4rhosigmalapl2,   0, dim->v4rhosigmalapl2  *np*sizeof(double));
      *v4rhosigmalapltau = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigmalapltau);
      libxc_memset(*v4rhosigmalapltau, 0, dim->v4rhosigmalapltau*np*sizeof(double));
      *v4rhosigmatau2    = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigmatau2);
      libxc_memset(*v4rhosigmatau2, 0, dim->v4rhosigmatau2*np*sizeof(double));
      *v4rholapl3        = (double *) libxc_malloc(sizeof(double)*np*dim->v4rholapl3);
      libxc_memset(*v4rholapl3,        0, dim->v4rholapl3       *np*sizeof(double));
      *v4rholapl2tau     = (double *) libxc_malloc(sizeof(double)*np*dim->v4rholapl2tau);
      libxc_memset(*v4rholapl2tau,     0, dim->v4rholapl2tau    *np*sizeof(double));
      *v4rholapltau2     = (double *) libxc_malloc(sizeof(double)*np*dim->v4rholapltau2);
      libxc_memset(*v4rholapltau2,     0, dim->v4rholapltau2    *np*sizeof(double));
      *v4rhotau3         = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhotau3);
      libxc_memset(*v4rhotau3,      0, dim->v4rhotau3     *np*sizeof(double));
      *v4sigma3lapl      = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma3lapl);
      libxc_memset(*v4sigma3lapl,      0, dim->v4sigma3lapl     *np*sizeof(double));
      *v4sigma3tau       = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma3tau);
      libxc_memset(*v4sigma3tau,    0, dim->v4sigma3tau   *np*sizeof(double));
      *v4sigma2lapl2     = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma2lapl2);
      libxc_memset(*v4sigma2lapl2,     0, dim->v4sigma2lapl2    *np*sizeof(double));
      *v4sigma2lapltau   = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma2lapltau);
      libxc_memset(*v4sigma2lapltau,   0, dim->v4sigma2lapltau  *np*sizeof(double));
      *v4sigma2tau2      = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma2tau2);
      libxc_memset(*v4sigma2tau2,   0, dim->v4sigma2tau2  *np*sizeof(double));
      *v4sigmalapl3      = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmalapl3);
      libxc_memset(*v4sigmalapl3,      0, dim->v4sigmalapl3     *np*sizeof(double));
      *v4sigmalapl2tau   = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmalapl2tau);
      libxc_memset(*v4sigmalapl2tau,   0, dim->v4sigmalapl2tau  *np*sizeof(double));
      *v4sigmalapltau2   = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmalapltau2);
      libxc_memset(*v4sigmalapl2tau,   0, dim->v4sigmalapl2tau  *np*sizeof(double));
      *v4sigmatau3       = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmatau3);
      libxc_memset(*v4sigmatau3,    0, dim->v4sigmatau3   *np*sizeof(double));
      *v4lapl4           = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapl4);
      libxc_memset(*v4lapl4,           0, dim->v4lapl4          *np*sizeof(double));
      *v4lapl3tau        = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapl3tau);
      libxc_memset(*v4lapl3tau,        0, dim->v4lapl3tau       *np*sizeof(double));
      *v4lapl2tau2       = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapl2tau2);
      libxc_memset(*v4lapl2tau2,       0, dim->v4lapl2tau2      *np*sizeof(double));
      *v4lapltau3        = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapltau3);
      libxc_memset(*v4lapltau3,        0, dim->v4lapltau3       *np*sizeof(double));
      *v4tau4            = (double *) libxc_malloc(sizeof(double)*np*dim->v4tau4);
      libxc_memset(*v4tau4,         0, dim->v4tau4        *np*sizeof(double));
    }
  }
#endif
#endif
#endif
#endif
}

void
xc_mgga_vars_free_all(double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{
  /* deallocate internal buffers */
  safe_free(zk);
#ifndef XC_DONT_COMPILE_VXC
  safe_free(vrho);
  safe_free(vsigma);
  safe_free(vlapl);
  safe_free(vtau);

#ifndef XC_DONT_COMPILE_FXC
  safe_free(v2rho2); safe_free(v2rhosigma); safe_free(v2rholapl); safe_free(v2rhotau);
  safe_free(v2sigma2); safe_free(v2sigmalapl); safe_free(v2sigmatau);
  safe_free(v2lapl2); safe_free(v2lapltau); safe_free(v2tau2);

#ifndef XC_DONT_COMPILE_KXC
  safe_free(v3rho3); safe_free(v3rho2sigma); safe_free(v3rho2lapl); safe_free(v3rho2tau);
  safe_free(v3rhosigma2); safe_free(v3rhosigmalapl); safe_free(v3rhosigmatau);
  safe_free(v3rholapl2); safe_free(v3rholapltau); safe_free(v3rhotau2);
  safe_free(v3sigma3); safe_free(v3sigma2lapl); safe_free(v3sigma2tau);
  safe_free(v3sigmalapl2); safe_free(v3sigmalapltau); safe_free(v3sigmatau2);
  safe_free(v3lapl3); safe_free(v3lapl2tau); safe_free(v3lapltau2); safe_free(v3tau3);

#ifndef XC_DONT_COMPILE_LXC
  safe_free(v4rho4); safe_free(v4rho3sigma); safe_free(v4rho3lapl); safe_free(v4rho3tau);
  safe_free(v4rho2sigma2); safe_free(v4rho2sigmalapl); safe_free(v4rho2sigmatau);
  safe_free(v4rho2lapl2); safe_free(v4rho2lapltau); safe_free(v4rho2tau2);
  safe_free(v4rhosigma3); safe_free(v4rhosigma2lapl); safe_free(v4rhosigma2tau);
  safe_free(v4rhosigmalapl2); safe_free(v4rhosigmalapltau); safe_free(v4rhosigmatau2);
  safe_free(v4rholapl3); safe_free(v4rholapl2tau); safe_free(v4rholapltau2); safe_free(v4rhotau3);
  safe_free(v4sigma4); safe_free(v4sigma3lapl); safe_free(v4sigma3tau); safe_free(v4sigma2lapl2);
  safe_free(v4sigma2lapltau); safe_free(v4sigma2tau2); safe_free(v4sigmalapl3); safe_free(v4sigmalapl2tau);
  safe_free(v4sigmalapltau2); safe_free(v4sigmatau3); safe_free(v4lapl4); safe_free(v4lapl3tau);
  safe_free(v4lapl2tau2); safe_free(v4lapltau3); safe_free(v4tau4);
#endif
#endif
#endif
#endif
}

void
xc_mgga_evaluate_functional(const xc_func_type *func, size_t np,
                            const double *rho, const double *sigma, const double *lapl, const double *tau,
                            double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{
  double *mzk = NULL;

  if(func->info->flags & XC_FLAGS_HAVE_EXC)
    mzk = zk;
  
  /* Evaluate the functional */
  switch(func->info->family){
  case XC_FAMILY_LDA:
  case XC_FAMILY_HYB_LDA:
    xc_lda(func, np, rho,
           mzk LDA_OUT_PARAMS_NO_EXC(XC_COMMA, ));
    break;
  case XC_FAMILY_GGA:
  case XC_FAMILY_HYB_GGA:
    xc_gga(func, np, rho, sigma,
           mzk GGA_OUT_PARAMS_NO_EXC(XC_COMMA, ));
    break;
  case XC_FAMILY_MGGA:
  case XC_FAMILY_HYB_MGGA:
    xc_mgga(func, np, rho, sigma, lapl, tau,
            mzk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ));
    break;
  }
}

/* initializes the mixing */
void
xc_deorbitalize_init(xc_func_type *p, int mgga_id, int ked_id)
{
  assert(p != NULL && p->func_aux == NULL);

  /* allocate structures needed for */
  p->n_func_aux = 2;
  p->func_aux   = (xc_func_type **) libxc_malloc(2*sizeof(xc_func_type *));

  p->func_aux[0] = (xc_func_type *) libxc_malloc(sizeof(xc_func_type));
  p->func_aux[1] = (xc_func_type *) libxc_malloc(sizeof(xc_func_type));

  xc_func_init (p->func_aux[0], mgga_id, p->nspin);
  xc_func_init (p->func_aux[1], ked_id,  p->nspin);
}

void
xc_deorbitalize_func_work(const xc_func_type *func, size_t np,
                     const double *rho, const double *sigma, const double *lapl, const double *tau,
                     double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{
  double *mrho, *msigma, *mlapl, *mtau;
  const double *null = NULL;
  double *mgga_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA *, mgga_);
  double *ked1_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA *, ked1_);
  double *ked2_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA *, ked2_);
  size_t ii;
  int order = -1;

  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;
  if(v3rho3 != NULL) order = 3;
  if(v4rho4 != NULL) order = 4;

  if(order < 0) return;

  /* prepare buffers that will hold the results from the individual functionals */
  mgga_zk MGGA_OUT_PARAMS_NO_EXC(=, mgga_ ) = NULL;
  ked1_zk MGGA_OUT_PARAMS_NO_EXC(=, ked1_ ) = NULL;
  ked2_zk MGGA_OUT_PARAMS_NO_EXC(=, ked2_ ) = NULL;

  /* allocate buffers */
  xc_mgga_vars_allocate_all(func->func_aux[0]->info->family, np, &(func->func_aux[0]->dim),
                       order >= 0, order >= 1, order >= 2, order >= 3, order >= 4,
                       &mgga_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, mgga_));
  xc_mgga_vars_allocate_all(func->func_aux[1]->info->family, np, &(func->func_aux[1]->dim),
                       order >= 0, order >= 1, order >= 2, order >= 3, order >= 4,
                       &ked1_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ked1_));

  if(func->nspin == XC_UNPOLARIZED){
    mtau   = (double *) libxc_malloc(sizeof(double)*np);
  }else{
    mrho   = (double *) libxc_malloc(2*sizeof(double)*np);
    msigma = (double *) libxc_malloc(3*sizeof(double)*np);
    mlapl  = (double *) libxc_malloc(2*sizeof(double)*np);
    mtau   = (double *) libxc_malloc(2*sizeof(double)*np);

    xc_mgga_vars_allocate_all(func->func_aux[1]->info->family, np, &(func->func_aux[1]->dim),
                         order >= 0, order >= 1, order >= 2, order >= 3, order >= 4,
                         &ked2_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ked2_));
  }

  /* evaluate the kinetic energy functional */
  if(func->nspin == XC_UNPOLARIZED){
    xc_mgga_evaluate_functional(func->func_aux[1], np, rho, sigma, lapl, NULL,
                           ked1_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ked1_));
  }else{
    for(ii=0; ii<np; ii++){
      mrho  [2*ii] = rho  [2*ii]; mrho  [2*ii+1] = 0.0;
      msigma[3*ii] = sigma[3*ii]; msigma[3*ii+1] = 0.0; msigma[3*ii+2] = 0.0;
      mlapl [2*ii] = lapl [2*ii]; mlapl [2*ii+1] = 0.0;
    }
    xc_mgga_evaluate_functional(func->func_aux[1], np, mrho, msigma, mlapl, NULL,
                           ked1_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ked1_));

    for(ii=0; ii<np; ii++){
      mrho  [2*ii] = rho  [2*ii + 1];
      msigma[3*ii] = sigma[3*ii + 2];
      mlapl [2*ii] = lapl [2*ii + 1];
    }
    xc_mgga_evaluate_functional(func->func_aux[1], np, mrho, msigma, mlapl, NULL,
                           ked2_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ked2_));

  }

  /* now evaluate the mgga functional */
  if(func->nspin == XC_UNPOLARIZED){
    for(ii=0; ii<np; ii++){
      mtau[ii] = rho[ii]*ked1_zk[ii];
    }
  }else{
    for(ii=0; ii<np; ii++){
      mtau[2*ii    ] = rho[2*ii    ]*ked1_zk[ii];
      mtau[2*ii + 1] = rho[2*ii + 1]*ked2_zk[ii];
    }
  }
  xc_mgga_evaluate_functional(func->func_aux[0], np, rho, sigma, lapl, mtau,
                         mgga_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, mgga_));

  /* now we have to combine the results */
  for(ii=0; ii<np; ii++){
    if(zk != NULL){
      *zk = *mgga_zk;
    }

#ifndef XC_DONT_COMPILE_VXC
    if(vrho != NULL){
#include "maple2c/deorbitalize_1.c"
    }
#ifndef XC_DONT_COMPILE_FXC
    if(v2rho2 != NULL){
#include "maple2c/deorbitalize_2.c"
    }
#ifndef XC_DONT_COMPILE_KXC
    if(v3rho3 != NULL){
#include "maple2c/deorbitalize_3.c"
    }
#ifndef XC_DONT_COMPILE_LXC
    if(v4rho4 != NULL){
#include "maple2c/deorbitalize_4.c"
    }
#endif
#endif
#endif
#endif

    internal_counters_mgga_next(&(func->dim), 0, &null, &null, &null, &null,
                                &zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ));
    internal_counters_mgga_next(&(func->func_aux[0]->dim), 0, &null, &null, &null, &null,
                                &mgga_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, mgga_));
    internal_counters_mgga_next(&(func->func_aux[1]->dim), 0, &null, &null, &null, &null,
                                &ked1_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ked1_));
    if(func->nspin == XC_POLARIZED){
      internal_counters_mgga_next(&(func->func_aux[1]->dim), 0, &null, &null, &null, &null,
                                  &ked2_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ked2_));
    }
  }

  /* move the counters back to zero and deallocate the memory */
  internal_counters_mgga_random(&(func->func_aux[0]->dim), -np, 0, &null, &null, &null, &null,
                                &mgga_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, mgga_));
  xc_mgga_vars_free_all(mgga_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, mgga_));

  internal_counters_mgga_random(&(func->func_aux[1]->dim), -np, 0, &null, &null, &null, &null,
                                &ked1_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ked1_));
  xc_mgga_vars_free_all(ked1_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ked1_));
  if(func->nspin == XC_POLARIZED){
    internal_counters_mgga_random(&(func->func_aux[1]->dim), -np, 0, &null, &null, &null, &null,
                                  &ked2_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ked2_));
    xc_mgga_vars_free_all(ked2_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ked2_));

    free(mrho); free(msigma); free(mlapl);
  }
  free(mtau);
}

static void
deorb_new(const xc_func_type *func, size_t np,
          const double *rho, const double *sigma, const double *lapl, const double *tau,
          xc_mgga_out_params *out)
{
  xc_deorbitalize_func_work(func, np, rho, sigma, lapl, tau,
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

xc_mgga_funcs_variants xc_deorbitalize_func =
  {
   {deorb_new, deorb_new, deorb_new, deorb_new, deorb_new},
   {deorb_new, deorb_new, deorb_new, deorb_new, deorb_new}
  };
                                               
