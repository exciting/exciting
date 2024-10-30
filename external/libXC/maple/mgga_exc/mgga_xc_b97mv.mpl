(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_xc_b97_mv_params *params;

  assert(p->params != NULL);
  params = (mgga_xc_b97_mv_params * )(p->params);
*)

b97mv_par_n := 5:

b97mv_gamma_x := 0.004:
b97mv_par_x := [
    [params_a_c_x[1], 0, 0],
    [params_a_c_x[2], 0, 1],
    [params_a_c_x[3], 0, 2],
    [params_a_c_x[4], 1, 0],
    [params_a_c_x[5], 1, 1]
]:

b97mv_gamma_ss := 0.2:
b97mv_par_ss := [
    [params_a_c_ss[1], 0, 0],
    [params_a_c_ss[2], 0, 2],
    [params_a_c_ss[3], 1, 0],
    [params_a_c_ss[4], 3, 2],
    [params_a_c_ss[5], 4, 2]
]:

b97mv_gamma_os := 0.006:
b97mv_par_os := [
    [params_a_c_os[1], 0, 0],
    [params_a_c_os[2], 0, 1],
    [params_a_c_os[3], 0, 3],
    [params_a_c_os[4], 1, 0],
    [params_a_c_os[5], 3, 2]
]:

$define lda_x_params
$include "lda_x.mpl"
$include "b97mv.mpl"

# This is the exchange term that also has to be added to the functional
b97mv_f_aux := (rs, z, xs0, xs1, ts0, ts1) ->
  + opz_pow_n( z,1)/2 * f_lda_x(rs*(2/(1 + z))^(1/3),  1)
    * b97mv_g(b97mv_gamma_x,  b97mv_wx_ss, b97mv_ux_ss, b97mv_par_x,  b97mv_par_n, xs0, ts0, 0)
  + opz_pow_n(-z,1)/2 * f_lda_x(rs*(2/(1 - z))^(1/3),  1)
    * b97mv_g(b97mv_gamma_x,  b97mv_wx_ss, b97mv_ux_ss, b97mv_par_x,  b97mv_par_n, xs1, ts1, 0):

f :=  (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  + b97mv_f_aux(rs, z, xs0, xs1, ts0, ts1)
  + b97mv_f(rs, z, xs0, xs1, ts0, ts1):
