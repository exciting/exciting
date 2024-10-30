(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* One should be able to simplify this by evaluating explicitly the Stoll
   decomposition of the exchange functional *)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

b97mv_ux    := (mgamma, x) -> mgamma*x^2/(1 + mgamma*x^2):

(* We define these here so that they can be redefined *)
b97mv_ux_ss := (mgamma, x) -> b97mv_ux(mgamma, x):
b97mv_ux_os := (mgamma, x) -> b97mv_ux(mgamma, x):

(* article uses t = 2 tau convention *)
b97mv_wx_ss := (t, dummy) -> (K_FACTOR_C - t)/(K_FACTOR_C + t):
b97mv_wx_os := (ts0, ts1) -> (K_FACTOR_C*(ts0 + ts1) - 2*ts0*ts1)/(K_FACTOR_C*(ts0 + ts1) + 2*ts0*ts1):

(* this is the generic Taylor expansion in wx and ux *)
b97mv_g := (mgamma, wx, ux, cc, n, xs, ts0, ts1) ->
  add(cc[i][1]*wx(ts0, ts1)^cc[i][2]*ux(mgamma, xs)^cc[i][3], i=1..n):

(* parallel contribution to correlation *)
b97mv_fpar  := (rs, z, xs0, xs1, ts0, ts1) ->
  + lda_stoll_par(f_pw, rs,  z,  1) * b97mv_g(b97mv_gamma_ss, b97mv_wx_ss, b97mv_ux_ss, b97mv_par_ss, b97mv_par_n, xs0, ts0, 0)
  + lda_stoll_par(f_pw, rs, -z, -1) * b97mv_g(b97mv_gamma_ss, b97mv_wx_ss, b97mv_ux_ss, b97mv_par_ss, b97mv_par_n, xs1, ts1, 0):

(* anti-parallel contribution *)
b97mv_fos := (rs, z, xs0, xs1, ts0, ts1) ->
  lda_stoll_perp(f_pw, rs, z)
  * b97mv_g(b97mv_gamma_os, b97mv_wx_os, b97mv_ux_os, b97mv_par_os, b97mv_par_n, sqrt(xs0^2 + xs1^2)/sqrt(2), ts0, ts1):

(* the functional is the sum of parallel and anti-parallel contributions *)
b97mv_f :=  (rs, z, xs0, xs1, ts0, ts1) ->
  + b97mv_fpar(rs, z, xs0, xs1, ts0, ts1)
  + b97mv_fos(rs, z, xs0, xs1, ts0, ts1):
