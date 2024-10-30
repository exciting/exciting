(*
 Copyright (C) 2020 Miguel Marques
               2024 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

(* eq 21 *)
f_ls := (rs, zeta) -> 0.0311*log(rs) - 0.048 + 0.009*rs*log(rs) - 0.01*rs:

(* eq 22: note that the paper appears to have a wrong overall sign;
this expression leads to a continuous interpolation and agreement with
Fig. 1 of the paper. It also appears that Fig. 1 has been scaled by 10 *)
f_ms := (rs, zeta) -> -0.06156 + 0.01898*log(rs):
(* eq 20 *)
f_hs := (rs, zeta) -> -0.438/rs + 1.325/rs^(3/2) - 1.47/rs^2 - 0.4/rs^(5/2):

f := (rs, zeta) -> my_piecewise5(
  rs < 0.7, f_ls(rs, zeta),  rs < 10,  f_ms(rs, zeta), f_hs(rs, zeta)):

