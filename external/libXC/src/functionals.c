/*
 Copyright (C) 2006-2021 M.A.L. Marques
               2015-2021 Susi Lehtola
               2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "xc.h"
#include "funcs_key.c"
#include <string.h>
#ifdef _MSC_VER
#define strcasecmp _stricmp
#define strncasecmp _strnicmp
#else
#include <strings.h>
#endif

extern xc_func_info_type
  *xc_lda_known_funct[],
  *xc_hyb_lda_known_funct[],
  *xc_gga_known_funct[],
  *xc_hyb_gga_known_funct[],
  *xc_mgga_known_funct[],
  *xc_hyb_mgga_known_funct[];


int xc_functional_get_number(const char *name)
{
  int ii;
  int key=-1;
  const char *p;

  /* Does name begin with xc_? */
  if(strncasecmp(name,"XC_",3) == 0) {
    p=name+3;
  } else {
    p=name;
  }

  for(ii=0;;ii++){
    if(xc_functional_keys[ii].number == -1)
      break;
    if(strcasecmp(xc_functional_keys[ii].name, p) == 0){
      key = xc_functional_keys[ii].number;
      break;
    }
  }

  return key;
}


char *xc_functional_get_name(int number)
{
  int ii;
  char *p;

  for(ii=0;;ii++){
    if(xc_functional_keys[ii].number == -1)
      return NULL;
    if(xc_functional_keys[ii].number == number) {
      /* return duplicated: caller has the responsibility to dealloc string.
         Do this the old way since strdup and strndup aren't C standard. */
      p = (char *) libxc_malloc(strlen(xc_functional_keys[ii].name) + 1);
      strcpy(p,xc_functional_keys[ii].name);
      return p;
    }
  }
}


int xc_family_from_id(int id, int *family, int *number)
{
  int ii;

  /* first let us check if it is an LDA */
  for(ii=0; xc_lda_known_funct[ii]!=NULL; ii++){
    if(xc_lda_known_funct[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_LDA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_LDA;
    }
  }

  for(ii=0; xc_hyb_lda_known_funct[ii]!=NULL; ii++){
    if(xc_hyb_lda_known_funct[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_HYB_LDA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_HYB_LDA;
    }
  }

  /* or is it a GGA? */
  for(ii=0; xc_gga_known_funct[ii]!=NULL; ii++){
    if(xc_gga_known_funct[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_GGA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_GGA;
    }
  }

  for(ii=0; xc_hyb_gga_known_funct[ii]!=NULL; ii++){
    if(xc_hyb_gga_known_funct[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_HYB_GGA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_HYB_GGA;
    }
  }

  /* or is it a meta GGA? */
  for(ii=0; xc_mgga_known_funct[ii]!=NULL; ii++){
    if(xc_mgga_known_funct[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_MGGA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_MGGA;
    }
  }

  for(ii=0; xc_hyb_mgga_known_funct[ii]!=NULL; ii++){
    if(xc_hyb_mgga_known_funct[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_HYB_MGGA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_HYB_MGGA;
    }
  }

  return XC_FAMILY_UNKNOWN;
}

int xc_number_of_functionals(void)
{
  int num;

  for(num=0;;num++){
    if(xc_functional_keys[num].number == -1)
      return num;
  }

  fprintf(stderr, "Internal error in functionals.c\n");
  exit(1);
}

int xc_maximum_name_length(void)
{
  int i, N, maxlen, tmp;

  N=xc_number_of_functionals();

  maxlen=0;
  for(i=0;i<N;i++){
    tmp=strlen(xc_functional_keys[i].name);
    if(tmp > maxlen) maxlen=tmp;
  }

  return maxlen;
}

static int compare_int(const void *a, const void *b) {
  return *(int *)a - *(int *) b;
}

void xc_available_functional_numbers(int *list)
{
  int ii, N;
  N=xc_number_of_functionals();
  for(ii=0;ii<N;ii++){
    list[ii]=xc_functional_keys[ii].number;
  }
  /* Sort list by id */
  qsort(list, N, sizeof(int), compare_int);
}

static int compare_func_names(const void *a, const void *b) {
  int ia, ib;
  int fama, famb;
  int hyba, hybb;

  ia = *(int *)a;
  ib = *(int *)b;

  /* First we sort by the family: LDAs, GGAs, meta-GGAs */
  fama = xc_family_from_id(xc_functional_keys[ia].number, NULL, NULL);
  famb = xc_family_from_id(xc_functional_keys[ib].number, NULL, NULL);
  if(fama < famb)
    return -1;
  else if(fama > famb)
    return 1;

  /* Then we sort by hybrid type: non-hybrids go first */
  hyba = (strncmp(xc_functional_keys[ia].name, "hyb_", 4) == 0);
  hybb = (strncmp(xc_functional_keys[ib].name, "hyb_", 4) == 0);
  if(!hyba && hybb)
    return -1;
  else if(hyba && !hybb)
    return 1;

  /* Last we sort by name */
  return strcmp(xc_functional_keys[ia].name, xc_functional_keys[ib].name);
}

void xc_available_functional_numbers_by_name(int *list)
{
  int ii, N;

  /* Arrange list of functional IDs by name */
  N=xc_number_of_functionals();
  for(ii=0;ii<N;ii++) {
    list[ii]=ii;
  }
  qsort(list, N, sizeof(int), compare_func_names);
  /* Map the internal list to functional IDs */
  for(ii=0;ii<N;ii++){
    list[ii]=xc_functional_keys[list[ii]].number;
  }
}

void xc_available_functional_names(char **list)
{
  int ii, N;
  int *idlist;

  /* Arrange list of functional IDs by name */
  N=xc_number_of_functionals();
  idlist=(int *) libxc_malloc(N*sizeof(int));
  for(ii=0;ii<N;ii++) {
    idlist[ii]=ii;
  }
  qsort(idlist, N, sizeof(int), compare_func_names);

  /* Now store the correctly ordered output */
  for(ii=0;ii<N;ii++) {
    strcpy(list[ii],xc_functional_keys[idlist[ii]].name);
  }

  /* Deallocate work array */
  libxc_free(idlist);
}

xc_func_type *xc_func_alloc(void)
{
  xc_func_type *func;

  func = (xc_func_type *) libxc_malloc (sizeof (xc_func_type));
  return func;
}

void xc_func_nullify(xc_func_type *func)
{
  assert(func != NULL);

  func->info       = NULL;
  func->nspin      = XC_UNPOLARIZED;

  func->n_func_aux = 0;
  func->func_aux   = NULL;
  func->mix_coef   = NULL;

  func->cam_omega  = 0.0;
  func->cam_alpha  = 0.0;
  func->cam_beta   = 0.0;

  func->nlc_b = func->nlc_C = 0.0;

  func->ext_params = NULL;
  func->params     = NULL;

  func->dens_threshold  = 0.0;
  func->zeta_threshold  = 0.0;
  func->sigma_threshold = 0.0;
  func->tau_threshold   = 0.0;
}

int xc_func_init(xc_func_type *func, int functional, int nspin)
{
  int number;

  assert(func != NULL);
  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);

  xc_func_nullify(func);

  /* initialize structure */
  func->nspin       = nspin;

  // we have to make a copy because the *_known_funct arrays live in
  // host memory (libxc_malloc instead returns memory than can be read
  // from GPU and CPU).
  xc_func_info_type * finfo = (xc_func_info_type *) libxc_malloc(sizeof(xc_func_info_type));

  // initialize the dimension structure
  libxc_memset(&(func->dim), 0, sizeof(xc_dimensions));
  switch(xc_family_from_id(functional, NULL, &number)){
  case(XC_FAMILY_LDA):
    *finfo = *xc_lda_known_funct[number];
    internal_counters_set_lda(func->nspin, &(func->dim));
    break;

  case(XC_FAMILY_HYB_LDA):
    *finfo = *xc_hyb_lda_known_funct[number];
    internal_counters_set_lda(func->nspin, &(func->dim));
    break;

  case(XC_FAMILY_GGA):
    *finfo = *xc_gga_known_funct[number];
    internal_counters_set_gga(func->nspin, &(func->dim));
    break;

  case(XC_FAMILY_HYB_GGA):
    *finfo = *xc_hyb_gga_known_funct[number];
    internal_counters_set_gga(func->nspin, &(func->dim));
    break;

  case(XC_FAMILY_MGGA):
    *finfo = *xc_mgga_known_funct[number];
    internal_counters_set_mgga(func->nspin, &(func->dim));
    break;

  case(XC_FAMILY_HYB_MGGA):
    *finfo = *xc_hyb_mgga_known_funct[number];
    internal_counters_set_mgga(func->nspin, &(func->dim));
    break;

  default:
    return -2; /* family not found */
  }
  func->info = finfo;

#ifdef XC_ENFORCE_FERMI_HOLE_CURVATURE
  /* Set the flag to enforce Fermi curvature by default */
  func->info->flags = func->info->flags | XC_FLAGS_ENFORCE_FHC;
#endif

  /* this is initialized for each functional from the info */
  func->dens_threshold = func->info->dens_threshold;
  /* the density and sigma cutoffs should be connected, especially in
     the case of kinetic energy functionals. This is the correct
     scaling */
  func->sigma_threshold = pow(func->info->dens_threshold, 4.0/3.0);

  /* these are reasonable defaults */
  func->zeta_threshold  = DBL_EPSILON;
  func->tau_threshold   = 1e-20;

  /* see if we need to initialize the functional */
  if(func->info->init != NULL)
    func->info->init(func);

  /* see if we need to initialize the external parameters */
  if(func->info->ext_params.n > 0) {
    func->ext_params = (double *) libxc_malloc(func->info->ext_params.n * sizeof(double));
    xc_func_set_ext_params(func, func->info->ext_params.values);

    /* sanity check external parameter names and descriptions */
    for(int i=0; i<func->info->ext_params.n; i++) {
      if(func->info->ext_params.names[i] == NULL) {
        char * name = xc_functional_get_name(functional);
        fprintf(stderr,"Internal error in %s: external parameter %i name is NULL\n",name,i);
        free(name);
        xc_func_end(func);
        return -1;
      }
      if(func->info->ext_params.descriptions[i] == NULL) {
        char * name = xc_functional_get_name(functional);
        fprintf(stderr,"Internal error in %s: external parameter %i description is NULL\n",name,i);
        free(name);
        xc_func_end(func);
        return -1;
      }
    }
  }

  return 0;
}


void xc_func_end(xc_func_type *func)
{
  assert(func != NULL && func->info != NULL);

  /* call internal termination routine */
  if(func->info->end != NULL)
    func->info->end(func);

  /* terminate any auxiliary functional */
  if(func->n_func_aux > 0){
    int ii;

    for(ii=0; ii<func->n_func_aux; ii++){
      xc_func_end(func->func_aux[ii]);
      libxc_free(func->func_aux[ii]);
    }
    libxc_free(func->func_aux);
  }

  /* deallocate coefficients for mixed functionals */
  if(func->mix_coef != NULL)
    libxc_free(func->mix_coef);

  /* deallocate any used parameter */
  if(func->ext_params != NULL)
    libxc_free(func->ext_params);
  if(func->params != NULL)
    libxc_free(func->params);

  libxc_free((void *) func->info);

  xc_func_nullify(func);
}

void  xc_func_free(xc_func_type *p)
{
  libxc_free(p);
}

const xc_func_info_type *xc_func_get_info(const xc_func_type *p)
{
  return p->info;
}

void xc_func_set_dens_threshold(xc_func_type *p, double t_dens)
{
  int ii;

  if(t_dens  > 0.0)
    p->dens_threshold = t_dens;

  for(ii=0; ii<p->n_func_aux; ii++) {
    xc_func_set_dens_threshold(p->func_aux[ii], t_dens);
  }
}

void xc_func_set_zeta_threshold(xc_func_type *p, double t_zeta)
{
  int ii;

  if(t_zeta  > 0.0)
    p->zeta_threshold = t_zeta;

  for(ii=0; ii<p->n_func_aux; ii++) {
    xc_func_set_zeta_threshold(p->func_aux[ii], t_zeta);
  }
}

void xc_func_set_sigma_threshold(xc_func_type *p, double t_sigma)
{
  int ii;

  if(t_sigma  > 0.0)
    p->sigma_threshold = t_sigma;

  for(ii=0; ii<p->n_func_aux; ii++) {
    xc_func_set_sigma_threshold(p->func_aux[ii], t_sigma);
  }
}

void xc_func_set_tau_threshold(xc_func_type *p, double t_tau)
{
  int ii;

  if(t_tau  > 0.0)
    p->tau_threshold = t_tau;

  for(ii=0; ii<p->n_func_aux; ii++) {
    xc_func_set_tau_threshold(p->func_aux[ii], t_tau);
  }
}

void  xc_func_set_fhc_enforcement(xc_func_type *p, int on)
{
  if(on)
    p->info->flags = p->info->flags | XC_FLAGS_ENFORCE_FHC;
  // Turn flag off if only it was on; otherwise ^ can flip it back on
  else if(p->info->flags & XC_FLAGS_ENFORCE_FHC)
    p->info->flags = p->info->flags ^ XC_FLAGS_ENFORCE_FHC;
}

void
xc_func_set_ext_params(xc_func_type *p, const double *ext_params)
{
  int ii;

  assert(p->info->ext_params.n > 0);
  /* We always store the current parameters in the functional instance */
  for(ii=0; ii<p->info->ext_params.n; ii++)
    p->ext_params[ii] = (ext_params[ii] == XC_EXT_PARAMS_DEFAULT) ? p->info->ext_params.values[ii] : ext_params[ii];
  /* Update the parameters */
  p->info->ext_params.set(p, p->ext_params);
}

void
xc_func_get_ext_params(const xc_func_type *p, double *ext_params)
{
  int ii;
  assert(p->info->ext_params.n > 0);
  assert(ext_params != NULL);
  for(ii=0; ii<p->info->ext_params.n; ii++)
    ext_params[ii] = p->ext_params[ii];
}

static int
xc_func_find_ext_params_name(const xc_func_type *p, const char *name) {
  int ii;
  assert(p != NULL && p->info->ext_params.n > 0);
  for(ii=0; ii<p->info->ext_params.n; ii++){
    if(strcmp(p->info->ext_params.names[ii], name) == 0) {
      return ii;
    }
  }
  return -1;
}

void
xc_func_set_ext_params_name(xc_func_type *p, const char *name, double par)
{
  int ii;

  ii = xc_func_find_ext_params_name(p, name);
  assert(ii>=0);
  p->ext_params[ii] = par;
  xc_func_set_ext_params(p, p->ext_params);
}

double
xc_func_get_ext_params_name(const xc_func_type *p, const char *name)
{
  int ii;

  ii = xc_func_find_ext_params_name(p, name);
  assert(ii>=0);
  return p->ext_params[ii];
}

double
xc_func_get_ext_params_value(const xc_func_type *p, int index)
{
  return p->ext_params[index];
}

/* returns the NLC parameters */
void xc_nlc_coef(const xc_func_type *p, double *nlc_b, double *nlc_C)
{
  assert(p!=NULL);

  *nlc_b = p->nlc_b;
  *nlc_C = p->nlc_C;
}
