(*
 Copyright (C) 2017 M.A.L. Marques
 Copyright (C) 2018 Susi Lehtola
 Copyright (C) 2024 Dogukan Yilmaz

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_msrpbel_params *params;

  assert(p->params != NULL);
  params = (mgga_x_msrpbel_params * ) (p->params);
*)

(* Equation (3) in the 2019 paper by Smeets et al. is missing the third power in the numerator *)
msrpbel_fa := a -> (1 - a^2)^3 / (1 + a^3 + params_a_b*a^6):
msrpbel_f0 := (p, c) -> 1 + params_a_kappa*(1 - exp(-(MU_GE*p + c)/params_a_kappa)):

msrpbel_alpha := (t,x) -> (t - x^2/8)/(K_FACTOR_C + params_a_eta*x^2/8):

msrpbel_f := (x, u, t) -> msrpbel_f0(X2S^2*x^2, 0) + \
  msrpbel_fa(msrpbel_alpha(t,x))*(msrpbel_f0(X2S^2*x^2, params_a_c) - msrpbel_f0(X2S^2*x^2, 0)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(msrpbel_f, rs, z, xs0, xs1, u0, u1, t0, t1):
