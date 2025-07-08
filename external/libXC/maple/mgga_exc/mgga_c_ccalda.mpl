(*
 Copyright (C) 2022 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_c_ccalda_params *params;

  assert(p->params != NULL);
  params = (mgga_c_ccalda_params * ) (p->params);
*)

$include "mgga_c_cc.mpl"

(* equation 7; note that CCaLDA uses non-spin-polarized quantities
unlike many other meta-GGAs *)
ccalda_alpha := (z, xt, ts0, ts1) ->
  (t_total(z, ts0, ts1) - xt^2/8) / K_FACTOR_C * 2^(2/3):

(* equation 10 *)
ccalda_f_alpha := a -> (1 + params_a_c) * a / (1 + params_a_c * a):

(* functional is defined by equation 9, use short hand to avoid mistakes *)
f_ccalda0 := (f_alpha, f_cc, f_pw) -> f_alpha*f_cc + (1-f_alpha)*f_pw:
(* with real arguments *)
f_ccalda := (rs, z, xt, ts0, ts1) -> f_ccalda0(ccalda_f_alpha(ccalda_alpha(z, xt, ts0, ts1)), f_cc(rs, z, xt, ts0, ts1), f_pw(rs, z)):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  f_ccalda(rs, z, xt, ts0, ts1):
