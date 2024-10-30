/*
 Copyright (C) 2006-2018 M.A.L. Marques
 Copyright (C) 2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/**
 * @file work_mgga.c
 * @brief This file is to be included in MGGA functionals.
 */

#ifdef XC_DEBUG
#define __USE_GNU 1
#include <fenv.h>
#endif

/* macro to simpligy accessing the variables */
#define VAR(var, ip, index)         var[ip*p->dim.var + index]
#define WORK_MGGA_(order, spin)     work_mgga_ ## order ## _ ## spin
#define WORK_MGGA_GPU_(order, spin) work_mgga_gpu_ ## order ## _ ## spin
#define FUNC_(order, spin)          func_     ## order ## _ ## spin

/* we need double escaping of the preprocessor macros */
#define WORK_MGGA(order, spin)     WORK_MGGA_(order, spin)
#define WORK_MGGA_GPU(order, spin) WORK_MGGA_GPU_(order, spin)
#define FUNC(order, spin)          FUNC_(order, spin)

#ifndef HAVE_CUDA

static void
WORK_MGGA(ORDER_TXT, SPIN_TXT)
(const XC(func_type) *p, size_t np,
 const double *rho, const double *sigma, const double *lapl, const double *tau,
 xc_mgga_out_params *out)
{

#ifdef XC_DEBUG
  /* This throws an exception when floating point errors are encountered */
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  size_t ip;
  double dens;
  double my_rho[2]={0.0, 0.0};
  double my_sigma[3]={0.0, 0.0, 0.0};
  double my_tau[2]={0.0, 0.0};

  for(ip = 0; ip < np; ip++){
    /* Screen low density */
    dens = (p->nspin == XC_POLARIZED) ? VAR(rho, ip, 0) + VAR(rho, ip, 1) : VAR(rho, ip, 0);
    if(dens < p->dens_threshold)
      continue;

    /* sanity check of input parameters */
    my_rho[0] = m_max(p->dens_threshold, VAR(rho, ip, 0));
    my_sigma[0] = m_max(p->sigma_threshold * p->sigma_threshold, VAR(sigma, ip, 0));

    /* Many functionals shamelessly divide by tau, so we set a reasonable threshold */
    /* skip all checks on tau for the kinetic functionals */
    if(p->info->flags & XC_FLAGS_NEEDS_TAU){
      my_tau[0] = m_max(p->tau_threshold, VAR(tau, ip, 0));
      if(p->info->flags & XC_FLAGS_ENFORCE_FHC) {
        /* The Fermi hole curvature 1 - xs^2/(8*ts) must be positive */
        my_sigma[0] = m_min(my_sigma[0], 8.0*my_rho[0]*my_tau[0]);
      }
    }
    /* lapl can have any values */

    if(p->nspin == XC_POLARIZED){
      double s_ave;

      my_rho[1] = m_max(p->dens_threshold, VAR(rho, ip, 1));
      my_sigma[2] = m_max(p->sigma_threshold * p->sigma_threshold, VAR(sigma, ip, 2));

      if(p->info->flags & XC_FLAGS_NEEDS_TAU){
        my_tau[1] = m_max(p->tau_threshold, VAR(tau, ip, 1));
        if(p->info->flags & XC_FLAGS_ENFORCE_FHC) {
          /* The Fermi hole curvature 1 - xs^2/(8*ts) must be positive */
          my_sigma[2] = m_min(my_sigma[2], 8.0*my_rho[1]*my_tau[1]);
        }
      }
      
      my_sigma[1] = VAR(sigma, ip, 1);
      s_ave = 0.5*(my_sigma[0] + my_sigma[2]);
      /* | grad n |^2 = |grad n_up + grad n_down|^2 > 0 */
      my_sigma[1] = (my_sigma[1] >= -s_ave ? my_sigma[1] : -s_ave);
      /* Since |grad n_up - grad n_down|^2 > 0 we also have */
      my_sigma[1] = (my_sigma[1] <= +s_ave ? my_sigma[1] : +s_ave);
    }

    FUNC(ORDER_TXT, SPIN_TXT)(p, ip, my_rho, my_sigma, &VAR(lapl, ip, 0), my_tau, out);

    /* check for NaNs */
#ifdef XC_DEBUG
    {
      const xc_dimensions *dim = &(p->dim);
      int ii, is_OK = 1;

      if(out->zk != NULL)
        is_OK = is_OK & isfinite(out->VAR(zk, ip, 0));

      if(out->vrho != NULL){
        for(ii=0; ii < dim->vrho; ii++)
          is_OK = is_OK && isfinite(out->VAR(vrho, ip, ii));
        for(ii=0; ii < dim->vsigma; ii++)
          is_OK = is_OK && isfinite(out->VAR(vsigma, ip, ii));
        if(p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)
          for(ii=0; ii < dim->vlapl; ii++)
            is_OK = is_OK && isfinite(out->VAR(vlapl, ip, ii));
        if(p->info->flags & XC_FLAGS_NEEDS_TAU)
          for(ii=0; ii < dim->vtau; ii++)
            is_OK = is_OK && isfinite(out->VAR(vtau, ip, ii));
      }

      if(!is_OK){
        printf("Problem in the evaluation of the functional\n");
        if(p->nspin == XC_UNPOLARIZED){
          printf("./xc-get_data %d 1 ", p->info->number);
          if(p->info->flags & (XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU))
            printf("%le 0.0 %le 0.0 0.0 %le 0.0 %le 0.0\n",
                   VAR(rho, ip, 0), VAR(sigma, ip, 0), VAR(lapl, ip, 0), VAR(tau, ip, 0));
          else if(p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)
            printf("%le 0.0 %le 0.0 0.0 %le 0.0 0.0 0.0\n",
                   VAR(rho, ip, 0), VAR(sigma, ip, 0), VAR(lapl, ip, 0));
          else
            printf("%le 0.0 %le 0.0 0.0 0.0 0.0 %le 0.0\n",
                   VAR(rho, ip, 0), VAR(sigma, ip, 0), VAR(tau, ip, 0));
        }else{
          printf("./xc-get_data %d 2 ", p->info->number);
          if(p->info->flags & (XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU))
            printf("%le %le %le %le %le %le %le %le %le\n",
                   VAR(rho, ip, 0), VAR(rho, ip, 1),
                   VAR(sigma, ip, 0), VAR(sigma, ip, 1), VAR(sigma, ip, 2),
                   VAR(lapl, ip, 0), VAR(lapl, ip, 1),
                   VAR(tau, ip, 0), VAR(tau, ip, 1));
          else if(p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)
            printf("%le %le %le %le %le %le %le 0.0 0.0\n",
                   VAR(rho, ip, 0), VAR(rho, ip, 1),
                   VAR(sigma, ip, 0), VAR(sigma, ip, 1), VAR(sigma, ip, 2),
                   VAR(lapl, ip, 0), VAR(lapl, ip, 1));
          else
            printf("%le %le %le %le %le 0.0 0.0 %le %le\n",
                   VAR(rho, ip, 0), VAR(rho, ip, 1),
                   VAR(sigma, ip, 0), VAR(sigma, ip, 1), VAR(sigma, ip, 2),
                   VAR(tau, ip, 0), VAR(tau, ip, 1));
        }
      }
    }
#endif
  }   /* for(ip) */

}

#else

__global__ static void
WORK_MGGA_GPU(ORDER_TXT, SPIN_TXT)
(const XC(func_type) *p, size_t np,
 const double *rho, const double *sigma, const double *lapl, const double *tau,
 xc_mgga_out_params *out)
{

  size_t ip = blockIdx.x * blockDim.x + threadIdx.x;
  double my_rho[2] = {0.0, 0.0};
  double my_sigma[3] = {0.0, 0.0, 0.0};
  double my_tau[2] = {0.0, 0.0};
  double dens;

  if(ip >= np) return;

  /* Screen small densities */
  dens = (p->nspin == XC_POLARIZED) ? VAR(rho, ip, 0) + VAR(rho, ip, 1) : VAR(rho, ip, 0);
  if(dens < p->dens_threshold)
    return;

  /* sanity check of input parameters */
  my_rho[0] = m_max(p->dens_threshold, VAR(rho, ip, 0));
  my_sigma[0] = m_max(p->sigma_threshold * p->sigma_threshold, VAR(sigma, ip, 0));
  
  /* Many functionals shamelessly divide by tau, so we set a reasonable threshold */
  if(p->info->flags & XC_FLAGS_NEEDS_TAU){
    my_tau[0] = m_max(p->tau_threshold, VAR(tau, ip, 0));
    if(p->info->flags & XC_FLAGS_ENFORCE_FHC) {
      /* The Fermi hole curvature 1 - xs^2/(8*ts) must be positive */
      my_sigma[0] = m_min(my_sigma[0], 8.0*my_rho[0]*my_tau[0]);
    }
  }
  /* lapl can have any values */
  if(p->nspin == XC_POLARIZED){
    double s_ave;

    my_rho[1]   = m_max(p->dens_threshold, VAR(rho, ip, 1));
    my_sigma[2] = m_max(p->sigma_threshold * p->sigma_threshold, VAR(sigma, ip, 2));
    if(p->info->flags & XC_FLAGS_NEEDS_TAU){
      my_tau[1] = m_max(p->tau_threshold, VAR(tau, ip, 1));
      if(p->info->flags & XC_FLAGS_ENFORCE_FHC) {
        /* The Fermi hole curvature 1 - xs^2/(8*ts) must be positive */
        my_sigma[2] = m_min(my_sigma[2], 8.0*my_rho[1]*my_tau[1]);
      }
    }
    my_sigma[1] = VAR(sigma, ip, 1);
    s_ave = 0.5*(my_sigma[0] + my_sigma[2]);
    /* | grad n |^2 = |grad n_up + grad n_down|^2 > 0 */
    my_sigma[1] = (my_sigma[1] >= -s_ave ? my_sigma[1] : -s_ave);
    /* Since |grad n_up - grad n_down|^2 > 0 we also have */
    my_sigma[1] = (my_sigma[1] <= +s_ave ? my_sigma[1] : +s_ave);
  }

  FUNC(ORDER_TXT, SPIN_TXT)(p, ip, my_rho, my_sigma, &VAR(lapl, ip, 0), my_tau, out);
}


static void
WORK_MGGA(ORDER_TXT, SPIN_TXT)
(const XC(func_type) *p, size_t np,
 const double *rho, const double *sigma, const double *lapl, const double *tau,
 xc_mgga_out_params *out)
{
  //make a copy of 'p' and 'out' since they might be in host-only memory
  XC(func_type) *pcuda = (XC(func_type) *) libxc_malloc(sizeof(XC(func_type)));
  xc_mgga_out_params *outcuda = (xc_mgga_out_params *) libxc_malloc(sizeof(xc_mgga_out_params));
  
  cudaMemcpy(pcuda, p, sizeof(XC(func_type)), cudaMemcpyHostToDevice);
  cudaMemcpy(outcuda, out, sizeof(xc_mgga_out_params), cudaMemcpyHostToDevice);

  size_t nblocks = np/CUDA_BLOCK_SIZE;
  if(np != nblocks*CUDA_BLOCK_SIZE) nblocks++;

  WORK_MGGA_GPU(ORDER_TXT, SPIN_TXT)<<<nblocks, CUDA_BLOCK_SIZE>>>
    (pcuda, np, rho, sigma, lapl, tau, outcuda);

  libxc_free(pcuda);
  libxc_free(outcuda);
}

#endif
