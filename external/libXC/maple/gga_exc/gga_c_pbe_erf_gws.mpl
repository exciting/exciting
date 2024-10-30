(*
 2018 Authored by Andrea Kreppel
 2022 Edited by Henryk Laqua

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 Short-range PBE correlation functional Goll/Werner/Stoll
 Goll, Werner, Stoll Phys. Chem. Chem. Phys. 7, (2005) 3917.
*)

(* type: gga_exc *)

(* prefix:
  gga_c_pbe_erf_gws_params *params;

  assert(p->params != NULL);
  params = (gga_c_pbe_erf_gws_params * )(p->params);
*)

$define lda_c_pw_params
$include "lda_c_pw_erf.mpl"
$include "lda_c_pw.mpl"

(*default parameters*)
$ifdef gga_c_pbe_erf_gws_params
params_a_beta  := 0.06672455060314922:
params_a_gamma := 0.031090690869654895034:
params_a_a_c    := 2.78
$endif

(*params*)
pbe_c_erf_gws_gamma := params_a_gamma:
pbe_c_erf_gws_beta_orig := params_a_beta:
pbe_c_erf_gws_a_c := params_a_a_c:

pbe_c_erf_gws_kS := (rs) -> (3/(4*Pi*n_total(rs)))^(1/3):

(*eq. (6)*)
pbe_c_erf_gws_beta := (rs, z) -> pbe_c_erf_gws_beta_orig * (lda_c_pw_erf_f(rs,z)/f_pw(rs,z))^pbe_c_erf_gws_a_c:

(* third eq. of eq. (6)*)
pbe_c_erf_gws_A := (rs, z) -> pbe_c_erf_gws_beta(rs, z)/(pbe_c_erf_gws_gamma*(exp(-lda_c_pw_erf_f(rs,z)/((mphi(z)^3)*pbe_c_erf_gws_gamma))-1)):

(* second eq. of eq. (6)*)
pbe_c_erf_gws_H := (rs, z, t) -> pbe_c_erf_gws_gamma*(mphi(z)^3)*ln(1 + pbe_c_erf_gws_beta(rs,z)*t^2/pbe_c_erf_gws_gamma*((1+pbe_c_erf_gws_A(rs,z)*t^2)/(1+pbe_c_erf_gws_A(rs,z)*t^2+pbe_c_erf_gws_A(rs,z)^2*t^4))):

(* first eq. of eq. (6)*)
f := (rs, z, xt, xs0, xs1) -> lda_c_pw_erf_f(rs,z) + pbe_c_erf_gws_H(rs, z, tt(rs,z,xt)):

