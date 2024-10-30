(*
 Copyright (C) 2024 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_bkl_params *params;

  assert(p->params != NULL);
  params = (gga_x_bkl_params * )(p->params);
*)

bkl_f0 := s -> 1 + params_a_gamma*params_a_kappa*(exp(-params_a_alpha*params_a_mu1*s^2) - exp(-params_a_beta*params_a_mu1*s^2)):
bkl_f  := x -> bkl_f0(X2S*x):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(bkl_f, rs, z, xs0, xs1):
