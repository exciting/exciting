/*
 2018 Authored by Andrea Kreppel
 2022 Edited by Henryk Laqua

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 Short-range (erfc) Perdew Wang correlation functional according to
 S. Paziani, S. Moroni, P. Gori-Giorgi, and G. B. Bachelet.,  Phys. Rev. B 73, 155111 (2006).
 DOI:10.1103/PhysRevB.73.155111
*/

#include "util.h"

#define  XC_LDA_C_PW_ERF                       654 /* Short ranged LDA correlation (erfc) */

static void
xc_lda_c_pw_erf_init(xc_func_type *p)
{
  xc_hyb_init_hybrid(p, 0.0);
}

static const char  *omega_names[]  = {"short_range_omega"};
static const char  *omega_desc[]   = {"short-range screening parameter"};
static const double omega_values[] = {0.5};

#include "maple2c/lda_exc/lda_c_pw_erf.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_pw_erf = {
  XC_LDA_C_PW_ERF,
  XC_CORRELATION,
  "Short ranged correlation LDA (erfc)",
  XC_FAMILY_LDA,
  {&xc_ref_Paziani2006_155111, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-13,
  {1, omega_names, omega_desc, omega_values, set_ext_params_cpy_omega},
  xc_lda_c_pw_erf_init, NULL,
  &work_lda, NULL, NULL
};
