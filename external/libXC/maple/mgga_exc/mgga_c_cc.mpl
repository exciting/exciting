(*
 Copyright (C) 2022 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

(* equation 2 in Lebeda 2022, tau_W / tau *)
cc_z := (z, xt, ts0, ts1) -> xt^2 / (8*t_total(z, ts0, ts1)):

(* equation 9 in Schmidt 2014, equation 6 in Lebeda 2022 *)
f_cc := (rs, z, xt, ts0, ts1) -> (1 - cc_z(z,xt,ts0,ts1)*z^2)*f_pw(rs, z):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  f_cc(rs, z, xt, ts0, ts1):
