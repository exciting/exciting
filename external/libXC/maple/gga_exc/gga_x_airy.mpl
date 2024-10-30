(*
 Copyright (C) 2017 M.A.L. Marques
               2024 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_airy_params *params;

  assert(p->params != NULL);
  params = (gga_x_airy_params * )(p->params);
*)

(* eq 18 *)
airy_f0 := s -> params_a_a1 * s^params_a_a2/(1 + params_a_a3 * s^params_a_a2)^params_a_a4 + (1 - params_a_a5*s^params_a_a6 + params_a_a7*s^params_a_a8)/(1 + params_a_a9*s^params_a_a10):

airy_f := x -> airy_f0(X2S*x):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(airy_f, rs, zeta, xs0, xs1):
