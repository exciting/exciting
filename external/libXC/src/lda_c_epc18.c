/*
 Copyright (C) 2023 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_C_EPC18_1     330  /* Electron-proton correlation 2018  */
#define XC_LDA_C_EPC18_2   331  /* Electron-proton correlation 2018 for proton affinities  */

typedef struct {
  double a, b, c;
} lda_c_epc18_params;

#define N_PAR 3
static const char  *names[N_PAR]  = {"_a", "_b", "_c"};
static const char  *desc[N_PAR]   = {"a", "b", "c"};

static const double par_epc18[N_PAR]   = {1.8, 0.1, 0.03};
static const double par_epc18_2[N_PAR] = {3.9, 0.5, 0.06};

static void
lda_c_epc18_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_c_epc18_params));
}

#include "maple2c/lda_exc/lda_c_epc18.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_epc18_1 = {
  XC_LDA_C_EPC18_1,
  XC_CORRELATION,
  "epc18-1: electron-proton correlation 2018",
  XC_FAMILY_LDA,
  {&xc_ref_Brorsen2018_044110, &xc_ref_epcnote, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_epc18, set_ext_params_cpy},
  lda_c_epc18_init, NULL,
  &work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_epc18_2 = {
  XC_LDA_C_EPC18_2,
  XC_CORRELATION,
  "epc18-2: electron-proton correlation 2018 for proton affinities",
  XC_FAMILY_LDA,
  {&xc_ref_Brorsen2018_044110, &xc_ref_epcnote, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_epc18_2, set_ext_params_cpy},
  lda_c_epc18_init, NULL,
  &work_lda, NULL, NULL
};
