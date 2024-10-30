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
  mgga_x_msb86bl_params *params;

  assert(p->params != NULL);
  params = (mgga_x_msb86bl_params * ) (p->params);
*)

(* Equation (3) in the 2019 paper by Smeets et al. is missing the third power in the numerator *)
msb86bl_fa := a -> (1 - a^2)^3 / (1 + a^3 + params_a_b*a^6):
msb86bl_f0 := (p, c) -> 1 + (MU_GE*p + c)/(1 + (MU_GE*p + c)/params_a_kappa)^(4/5):

msb86bl_alpha := (t,x) -> (t - x^2/8)/(K_FACTOR_C + params_a_eta*x^2/8):

msb86bl_f := (x, u, t) -> msb86bl_f0(X2S^2*x^2, 0) + \
  msb86bl_fa(msb86bl_alpha(t,x))*(msb86bl_f0(X2S^2*x^2, params_a_c) - msb86bl_f0(X2S^2*x^2, 0)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(msb86bl_f, rs, z, xs0, xs1, u0, u1, t0, t1):
