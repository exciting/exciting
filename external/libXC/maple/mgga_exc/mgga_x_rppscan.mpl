(*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_rppscan_params *params;

  assert(p->params != NULL);
  params = (mgga_x_rppscan_params * )(p->params);
*)

$include "mgga_x_scan.mpl"
$include "mgga_x_r2scan.mpl"
$include "mgga_x_rscan.mpl"

(* r++SCAN is obtained from rSCAN by replacing the definition of alpha *)
rscan_alpha := (rs, z, x, t) -> r2scan_alpha(x, t):
