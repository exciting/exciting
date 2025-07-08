(*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

(* prefix:
  mgga_x_ft98_params *params;

  assert(p->params != NULL);
  params = (mgga_x_ft98_params * ) (p->params);
*)

(* eq. 23 *)
ft98_f1 := xi -> (1 + params_a_a1*xi)^(1/2) / (1 + params_a_b1*xi)^(3/4):

(* eq. 24a *)
ft98_f2 := (xi, chi) -> (1 + params_a_a2*ft98_q1(xi, chi)) * (1 + ft98_q2(xi, chi)) / (1 + (2^(1/3) - 1) * ft98_q2(xi, chi))^3:
(* eq. 24b *)
ft98_q1 := (xi, chi) -> (xi-chi)^2 / (1 + xi)^2:

(* eq. 24c *)
ft98_q2_orig := q3 -> 1/(q3 + (1 + q3^2)^(1/2)):

(* handle small-|q| and very negative q regions separately.

   The small-q expansion is
   1 - q + 1/2 q^2 - 1/8 q^4 + ...

   which is accurate to epsilon when q ~ epsilon^(1/4)
*)
ft98_q2term_smallq := q3 -> eval(convert(taylor(ft98_q2_orig(qt), qt = 0, 9), polynom), qt=q3):
ft98_q2_cutoff_smallq := DBL_EPSILON^(1/4):

(* when q->-infty, there's danger for overflow.

   The expansion at -infty is

   -2q - 1/(2q) + 1/(8q^3) - ...

   which will be accurate when q < -epsilon^(-1/4)
*)
ft98_q2term_minfty := q3 -> eval(convert(taylor(ft98_q2_orig(qt), qt = -infinity, 9), polynom), qt=q3):
ft98_q2_cutoff_minfty := -DBL_EPSILON^(-1/4):

(* assemble the result. Linting the argument for the original argument
   is tricky since the arguments have several branches:

   q <= minfty: ft98_q2term_minfty
   minfty <= q <= -smallq: ft98_q2_orig
   -smallq <= q <= +smallq: ft98_q2_smallq
   smallq <= q: ft98_q2_orig
*)
ft98_q20 := q3 -> my_piecewise5(
      q3 < ft98_q2_cutoff_minfty, ft98_q2term_minfty(q3),
      m_abs(q3) < ft98_q2_cutoff_smallq, ft98_q2term_smallq(q3),
      ft98_q2_orig(m_max(q3, ft98_q2_cutoff_minfty))
  ):
ft98_q2 := (xi, chi) -> ((params_a_b2^2+1)^(1/2) - params_a_b2)*ft98_q20(ft98_q3(xi, chi)):

(* eq. 24d *)
ft98_q3 := (xi, chi) -> xi^2 - chi^2 - params_a_b2:

(* eq. 12 *)
ft98_f0 := (xi, chi) ->
  sqrt((1 + params_a_a*ft98_f1(xi)*xi + params_a_b*ft98_f2(xi, chi)*(xi-chi)^2)/(1 + 36*LDA_X_FACTOR^2*params_a_b*xi)):

(*
   eq. 4: the xi variable is simply libxc's x^2
   eq. 5: the chi variable is libxc's u variable
*)
ft98_f := (x, u, t) -> ft98_f0(x^2, u):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(ft98_f, rs, z, xs0, xs1, u0, u1, t0, t1):
