(*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_r4scan_params *params;

  assert(p->params != NULL);
  params = (mgga_x_r4scan_params * )(p->params);
*)

$include "mgga_x_r2scan.mpl"

(* r4SCAN is obtained from r2SCAN by replacing the enhancement factor *)

(* df2 ~= -0.9353000875519996 *)
df2 := ff -> add(i*ff[9-i], i=1..8):
(* df4 ~=  0.8500359204920018 *)
df4 := ff -> add((i-1)*(i-2)*ff[9-i], i=2..8):

(* eq 53 *)
Cn := (20/27 + params_a_eta*5/3):
(* eq 61 *)
Caa := ff -> 73/5000 - df4(ff)/2*(scan_h0x-1):
(* eq 62 *)
Cpa := ff -> 511/13500 - 73/1500*params_a_eta - df2(ff)*(Cn*C2(ff)+MU_GE):
(* eq 63 *)
Cpp := ff -> 146/2025*(params_a_eta*3/4 + 2/3)^2 - 73/405*(params_a_eta*3/4 + 2/3) + (Cn * C2(ff) + MU_GE)^2 / params_a_k1:

(* eq 59 *)
r4scan_dF := (ff, p, a) -> (C2(ff) * ((1-a)-Cn*p) + Caa(ff)*(1-a)^2 + Cpa(ff)*p*(1-a) + Cpp(ff)*p^2)*r4scan_dFdamp(p,a):
(* eq 60 *)
r4scan_dFdamp := (p, a) -> 2*a^2/(1+a^4) * exp(-(1-a)^2/params_a_da4^2 - p^2/params_a_dp4^4):

r4scan_f := (x, u, t) -> (scan_h1x(r2scan_x(scan_p(x), rscan_fx)) + r2scan_f_alpha(r2scan_alpha(x, t), rscan_fx) * (scan_h0x - scan_h1x(r2scan_x(scan_p(x), rscan_fx))) + r4scan_dF(rscan_fx, scan_p(x), r2scan_alpha(x, t)))*scan_gx(x):
f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(r4scan_f, rs, z, xs0, xs1, u0, u1, t0, t1):
