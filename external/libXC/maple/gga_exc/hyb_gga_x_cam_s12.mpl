(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  hyb_gga_x_cam_s12_params *params;

  assert(p->params != NULL);
  params = (hyb_gga_x_cam_s12_params * )(p->params);
*)

$include "gga_x_s12.mpl"
$include "gga_x_ityh.mpl"

(* we remove the params_a_bx parameter from the definition of s12g_f *)
ityh_enhancement := xs  -> s12g_f(xs)/params_a_bx:

cam_s12_f := (rs, z, xs) -> ityh_enhancement(xs) *
  (1 - p_a_cam_alpha - p_a_cam_beta*ityh_f_aa(rs, z, xs)):

f := (rs, z, xt, xs0, xs1) -> gga_exchange_nsp(cam_s12_f, rs, z, xs0, xs1):
