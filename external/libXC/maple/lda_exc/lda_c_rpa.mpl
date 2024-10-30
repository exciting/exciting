(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

(* Carr & Maradudin 1964, eq 27, dividing by two to convert from Rydberg to Hartree *)
f := (rs, zeta) -> ((0.0622*log(rs) - 0.096) + rs*(0.018*log(rs) - 0.036))/2:
