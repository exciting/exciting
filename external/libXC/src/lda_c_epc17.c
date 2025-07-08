/*
 Copyright (C) 2023 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_C_EPC17     328  /* Electron-proton correlation 2017 a.k.a. epc17-1   */
#define XC_LDA_C_EPC17_2   329  /* Electron-proton correlation 2017 for proton affinities  */

typedef struct {
  double a, b, c;
} lda_c_epc17_params;

#define N_PAR 3
static const char  *names[N_PAR]  = {"_a", "_b", "_c"};
static const char  *desc[N_PAR]   = {"a", "b", "c"};

static const double par_epc17[N_PAR]   = {2.35, 2.4, 3.2};
static const double par_epc17_2[N_PAR] = {2.35, 2.4, 6.6};

static void
lda_c_epc17_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_c_epc17_params));
}

#include "maple2c/lda_exc/lda_c_epc17.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_epc17 = {
  XC_LDA_C_EPC17,
  XC_CORRELATION,
  "epc17(-1): electron-proton correlation 2017",
  XC_FAMILY_LDA,
  {&xc_ref_Yang2017_114113, &xc_ref_epcnote, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_epc17, set_ext_params_cpy},
  lda_c_epc17_init, NULL,
  &work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_epc17_2 = {
  XC_LDA_C_EPC17_2,
  XC_CORRELATION,
  "epc17-2: electron-proton correlation 2017 for proton affinities",
  XC_FAMILY_LDA,
  {&xc_ref_Brorsen2017_3488, &xc_ref_epcnote, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_epc17_2, set_ext_params_cpy},
  lda_c_epc17_init, NULL,
  &work_lda, NULL, NULL
};
