(*
 Copyright (C) 2017 M.A.L. Marques
               2024 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_hcth_a_params *params;

  assert(p->params != NULL);
  params = (gga_x_hcth_a_params * )(p->params);
*)

(* eq 31 *)
hcth_b88x := (beta,x) -> beta*x^2/(1+params_a_gamma*beta*x*arcsinh(x)):
(* eq 30 *)
hcth_a_f := x -> params_a_c0 - params_a_c1/X_FACTOR_C*hcth_b88x(params_a_beta,x) - params_a_c2/X_FACTOR_C*eval(diff(hcth_b88x(beta,x),beta),beta=params_a_beta):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(hcth_a_f, rs, z, xs0, xs1):
