(*
 Copyright (C) 2017 M.A.L. Marques
               2022 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$include "gga_c_tca.mpl"

(* 4th order Taylor expansion of sinc(x) at x=0 *)
sinc := x -> sin(x)/x:
sinc_taylor := x -> eval(convert(taylor(sinc(y),y=0,9),polynom),y=x):

(* Switch to Taylor expansion when x^4 = epsilon. Remember that we
need further orders in the Taylor series to make the derivatives
accurate as well. *)
sinc_cutoff := DBL_EPSILON^(1/4):
msinc := x -> my_piecewise3(x <= sinc_cutoff, sinc_taylor(x), sinc(m_max(x, sinc_cutoff))):

revtca_aa := Pi*(9*Pi/4)^(1/3):
revtca_fD := (rs, z, s) -> 1 - z^4*(1 - msinc(revtca_aa*s/rs)^2):

f := (rs, z, xt, xs0, xs1) ->
  f_tcs(rs, z, xt)*revtca_fD(rs, z, X2S*2^(1/3)*xt):
