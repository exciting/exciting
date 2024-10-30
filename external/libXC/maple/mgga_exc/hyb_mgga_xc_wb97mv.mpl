(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  hyb_mgga_xc_wb97_mv_params *params;

  assert(p->params != NULL);
  params = (hyb_mgga_xc_wb97_mv_params * )(p->params);
*)

b97mv_par_n := 6:

b97mv_gamma_x := 0.004:
b97mv_par_x := [
    [  params_a_c_x[1], 0, 0],
    [  params_a_c_x[2], 0, 1],
    [  params_a_c_x[3], 1, 0],
    [  0,   0, 0],
    [  0,   0, 0],
    [  0,   0, 0]
]:

b97mv_gamma_ss := 0.2:
b97mv_par_ss := [
    [  params_a_c_ss[1],  0, 0],
    [  params_a_c_ss[2],  0, 4],
    [  params_a_c_ss[3],  1, 0],
    [  params_a_c_ss[4],  2, 0],
    [  params_a_c_ss[5],  4, 3],
    [  0,    0, 0]
]:

b97mv_gamma_os := 0.006:
b97mv_par_os := [
    [  params_a_c_os[1],  0, 0],
    [  params_a_c_os[2],  1, 0],
    [  params_a_c_os[3],  2, 0],
    [  params_a_c_os[4],  2, 1],
    [  params_a_c_os[5],  6, 0],
    [  params_a_c_os[6],  6, 1]
]:

$include "lda_x_erf.mpl"
$include "b97mv.mpl"

wb97mv_f := (rs, z, xs0, xs1, ts0, ts1) ->
   my_piecewise3(screen_dens_zeta(rs,  z), 0, (1 + z)/2 * lda_x_erf_spin(rs*(2/(1 + z))^(1/3),  1)
    * b97mv_g(b97mv_gamma_x,  b97mv_wx_ss, b97mv_ux_ss, b97mv_par_x,  b97mv_par_n, xs0, ts0, 0))
+  my_piecewise3(screen_dens_zeta(rs, -z), 0, (1 - z)/2 * lda_x_erf_spin(rs*(2/(1 - z))^(1/3),  1)
    * b97mv_g(b97mv_gamma_x,  b97mv_wx_ss, b97mv_ux_ss, b97mv_par_x,  b97mv_par_n, xs1, ts1, 0)):

f :=  (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  wb97mv_f(rs, z, xs0, xs1, ts0, ts1) +
  b97mv_f(rs, z, xs0, xs1, ts0, ts1):
