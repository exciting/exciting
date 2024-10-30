(*
 Copyright (C) 2023 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_c_epc18_params *params;

  assert(p->params != NULL);
  params = (lda_c_epc18_params * )(p->params);
*)

(* Equation (11) *)
epc18_beta := (rho_e, rho_p) -> rho_e^(1/3) + rho_p^(1/3):
(* Equation (7) *)
epc18_E := (rho_e, rho_p) -> - (rho_e*rho_p) / (params_a_a - params_a_b*epc18_beta(rho_e, rho_p)^3 + params_a_c * epc18_beta(rho_e, rho_p)^6):

(* Energy density using equation (11); require significant proton and electron density *)
f_epc18 := (rs, zeta) ->
        my_piecewise3(screen_dens(rs,zeta) and screen_dens(rs,-zeta), 0,
        epc18_E(n_spin(rs, z_thr(zeta)), n_spin(rs, z_thr(-zeta))) / n_total(rs)):

f := (rs, zeta) -> f_epc18(rs, zeta):
