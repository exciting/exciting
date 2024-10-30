/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_lda.c"
#include "funcs_hyb_lda.c"

void
xc_lda_sanity_check(const xc_func_info_type *info, int order, xc_lda_out_params *out)
{
  /* sanity check */
  if(order < 0 || order > 4){
    fprintf(stderr, "Order of derivatives '%d' not implemented\n",
	    order);
    exit(1);
  }

  if(out->zk != NULL && !(info->flags & XC_FLAGS_HAVE_EXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of Exc\n",
	    info->name);
    exit(1);
  }

  if(out->vrho != NULL && !(info->flags & XC_FLAGS_HAVE_VXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of vxc\n",
	    info->name);
    exit(1);
  }

  if(out->v2rho2 != NULL && !(info->flags & XC_FLAGS_HAVE_FXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of fxc\n",
	    info->name);
    exit(1);
  }

  if(out->v3rho3 != NULL && !(info->flags & XC_FLAGS_HAVE_KXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of kxc\n",
	    info->name);
    exit(1);
  }
}


void
xc_lda_initalize(const xc_func_type *func, size_t np, xc_lda_out_params *out)
{
  const xc_dimensions *dim = &(func->dim);

  /* initialize output */
  if(out->zk != NULL)
    libxc_memset(out->zk,     0, np*sizeof(double)*dim->zk);

  if(out->vrho != NULL)
    libxc_memset(out->vrho,   0, np*sizeof(double)*dim->vrho);

  if(out->v2rho2 != NULL)
    libxc_memset(out->v2rho2, 0, np*sizeof(double)*dim->v2rho2);

  if(out->v3rho3 != NULL)
    libxc_memset(out->v3rho3, 0, np*sizeof(double)*dim->v3rho3);

  if(out->v4rho4 != NULL)
    libxc_memset(out->v4rho4, 0, np*sizeof(double)*dim->v4rho4);
}


/* get the lda functional */
void
xc_lda_new(const xc_func_type *func, int order, size_t np, const double *rho,
       xc_lda_out_params *out)
{
  xc_lda_sanity_check(func->info, order, out);
  xc_lda_initalize(func, np, out);

  /* call the LDA routines */
  if(func->info->lda != NULL){
    if(func->nspin == XC_UNPOLARIZED){
      if(func->info->lda->unpol[order] != NULL)
        func->info->lda->unpol[order](func, np, rho, out);
    }else{
      if(func->info->lda->pol[order] != NULL)
        func->info->lda->pol[order](func, np, rho, out);
    }
  }

  if(func->mix_coef != NULL)
    xc_mix_func(func, np, rho, NULL, NULL, NULL, out->zk, out->vrho, NULL, NULL, NULL,
                out->v2rho2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                out->v3rho3, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                out->v4rho4, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL);
}

/* old API */
void
xc_lda(const xc_func_type *p, size_t np, const double *rho, double *zk, double *vrho, double *v2rho2, double *v3rho3, double *v4rho4)
{
  int order = -1;

  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;
  if(v3rho3 != NULL) order = 3;
  if(v4rho4 != NULL) order = 4;

  if(order < 0) return;

  xc_lda_out_params out;
  libxc_memset(&out, 0, sizeof(xc_lda_out_params));
  out.zk     = zk;
  out.vrho   = vrho;
  out.v2rho2 = v2rho2;
  out.v3rho3 = v3rho3;
  out.v4rho4 = v4rho4;

  xc_lda_new(p, order, np, rho, &out);
}


/* specializations */
void
xc_lda_exc(const xc_func_type *p, size_t np, const double *rho, double *zk)
{
  xc_lda_out_params out;
  libxc_memset(&out, 0, sizeof(xc_lda_out_params));
  out.zk   = zk;

  xc_lda_new(p, 0, np, rho, &out);
}

void
xc_lda_exc_vxc(const xc_func_type *p, size_t np, const double *rho, double *zk, double *vrho)
{
  xc_lda_out_params out;
  libxc_memset(&out, 0, sizeof(xc_lda_out_params));
  out.zk   = zk;
  out.vrho = vrho;

  xc_lda_new(p, 1, np, rho, &out);
}

void
xc_lda_exc_vxc_fxc(const xc_func_type *p, size_t np, const double *rho, double *zk, double *vrho, double *v2rho2)
{
  xc_lda_out_params out;
  libxc_memset(&out, 0, sizeof(xc_lda_out_params));
  out.zk     = zk;
  out.vrho   = vrho;
  out.v2rho2 = v2rho2;

  xc_lda_new(p, 2, np, rho, &out);
}

void
xc_lda_vxc_fxc(const xc_func_type *p, size_t np, const double *rho, double *vrho, double *v2rho2)
{
  xc_lda_out_params out;
  libxc_memset(&out, 0, sizeof(xc_lda_out_params));
  out.vrho   = vrho;
  out.v2rho2 = v2rho2;

  xc_lda_new(p, 2, np, rho, &out);
}

void
xc_lda_exc_vxc_fxc_kxc(const xc_func_type *p, size_t np, const double *rho, double *zk, double *vrho, double *v2rho2, double *v3rho3)
{
  xc_lda_out_params out;
  libxc_memset(&out, 0, sizeof(xc_lda_out_params));
  out.zk     = zk;
  out.vrho   = vrho;
  out.v2rho2 = v2rho2;
  out.v3rho3 = v3rho3;

  xc_lda_new(p, 3, np, rho, &out);
}

void
xc_lda_vxc_fxc_kxc(const xc_func_type *p, size_t np, const double *rho, double *vrho, double *v2rho2, double *v3rho3)
{
  xc_lda_out_params out;
  libxc_memset(&out, 0, sizeof(xc_lda_out_params));
  out.vrho   = vrho;
  out.v2rho2 = v2rho2;
  out.v3rho3 = v3rho3;

  xc_lda_new(p, 3, np, rho, &out);
}

void
xc_lda_vxc(const xc_func_type *p, size_t np, const double *rho, double *vrho)
{
  xc_lda_out_params out;
  libxc_memset(&out, 0, sizeof(xc_lda_out_params));
  out.vrho   = vrho;

  xc_lda_new(p, 1, np, rho, &out);
}

void
xc_lda_fxc(const xc_func_type *p, size_t np, const double *rho, double *v2rho2)
{
  xc_lda_out_params out;
  libxc_memset(&out, 0, sizeof(xc_lda_out_params));
  out.v2rho2 = v2rho2;

  xc_lda_new(p, 2, np, rho, &out);
}

void
xc_lda_kxc(const xc_func_type *p, size_t np, const double *rho, double *v3rho3)
{
  xc_lda_out_params out;
  libxc_memset(&out, 0, sizeof(xc_lda_out_params));
  out.v3rho3 = v3rho3;

  xc_lda_new(p, 3, np, rho, &out);
}

void
xc_lda_lxc(const xc_func_type *p, size_t np, const double *rho, double *v4rho4)
{
  xc_lda_out_params out;
  libxc_memset(&out, 0, sizeof(xc_lda_out_params));
  out.v4rho4 = v4rho4;

  xc_lda_new(p, 4, np, rho, &out);
}
