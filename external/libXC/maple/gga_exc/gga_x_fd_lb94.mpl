(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_fd_lb94_params *params;

  assert(p->params != NULL);
  params = (gga_x_fd_lb94_params * )(p->params);
*)

(* replace: "fd_int0\(" -> "xc_integrate(func0, &(params->beta), 0.0, " *)
(* replace: "fd_int1\(" -> "xc_integrate(func1, &(params->beta), 0.0, " *)

`diff/fd_int0` := proc(g, x) diff(g, x) * fd_f_inter(0, g) end proc:
`diff/fd_int1` := proc(g, x) diff(g, x) * fd_f_inter(1, g) end proc:

fd_beta := params_a_beta:
fd_csi  := 2^(1/3):

fd_f_inter := (n, t) -> -3/4 * fd_beta*fd_csi*log(t)^n /
  (1 + 3*fd_beta*fd_csi*t*log(fd_csi*t + sqrt((fd_csi*t)^2 + 1))):

fd_f := x -> 1 - x/X_FACTOR_C*(fd_int0(x/fd_csi)*log(x/fd_csi) - fd_int1(x/fd_csi)):

f  := (rs, z, xt, xs0, xs1) -> gga_exchange(fd_f, rs, z, xs0, xs1):
