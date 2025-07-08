(*
 Copyright (C) 2023 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_c_epc17_params *params;

  assert(p->params != NULL);
  params = (lda_c_epc17_params * )(p->params);
*)

(* Equation (13) *)
epc17_E := rho_ep -> -rho_ep / (params_a_a - params_a_b*sqrt(rho_ep) + params_a_c * rho_ep):

(* Energy density; require significant proton and electron density *)
f_epc17 := (rs, zeta) ->
        my_piecewise3(screen_dens(rs,zeta) and screen_dens(rs,-zeta), 0,
        epc17_E(n_spin(rs, z_thr(zeta)) * n_spin(rs, z_thr(-zeta))) / n_total(rs)):

f := (rs, zeta) -> f_epc17(rs, zeta):
