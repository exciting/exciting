(*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  hyb_gga_xc_case21_params *params;

  assert(p->params != NULL);
  params = (hyb_gga_xc_case21_params * )(p->params);
*)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"
$include "lda_x.mpl"

(* Teach maple to differentiate the B spline functions *)
`diff/xbspline` := (u, ider, params, x) -> xbspline(u, ider+1, params) * diff(u, x):
`diff/cbspline` := (u, ider, params, x) -> cbspline(u, ider+1, params) * diff(u, x):

(* text after eq 3, B86-type descriptor for exchange *)
case21_ux0 := s -> params_a_gammax*s^2/(1 + params_a_gammax*s^2):
(* enhancement function, eq 6 *)
case21_fx := x -> xbspline(case21_ux0(X2S*x), 0, params):
(* exchange energy, eq 3 *)
case21_Ex := (rs, z, xs0, xs1) -> gga_exchange(case21_fx, rs, z, xs0, xs1):

(* eq 5 *)
case21_t := (rs, z, xs0, xs1) -> (Pi/3)^(1/6) * (xs0*n_spin(rs, z)^(4/3) + xs1*n_spin(rs, -z)^(4/3))/( 4 * n_total(rs)^(7/6) * mphi(z)):
(* text before eq 5 *)
case21_uc := (rs, z, xs0, xs1) -> (-mphi(z)^3 * case21_t(rs, z, xs0, xs1)^2) / (-mphi(z)^3 * case21_t(rs, z, xs0, xs1)^2 + params_a_gammac*f_pw(rs, z)):
(* correlation energy, eqs 4 and 6 *)
case21_Ec := (rs, z, xs0, xs1) -> cbspline(case21_uc(rs, z, xs0, xs1), 0, params)*f_pw(rs, z):

(* whole functional *)
f_case21 := (rs, z, xt, xs0, xs1) ->
   (1-params_a_ax)*case21_Ex(rs, z, xs0, xs1) + case21_Ec(rs, z, xs0, xs1):

f  := (rs, z, xt, xs0, xs1) -> f_case21(rs, z, xt, xs0, xs1):
