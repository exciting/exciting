(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

(* replace: "mbrxc_x\(" -> "xc_mgga_x_mbrxc_get_x(" *)

mbrxc_a1 := 0.074746:
mbrxc_a2 := 0.147:
mbrxc_a3 := 0.0032:

(* This is the derivative of f = (1+x)^(5/3)*exp(-2/3*x)/(x - 3) = (32*Pi)^(2/3)/(6*Q) *)
mbrxc_aux_dfdx := x -> -2/3 * (1 + x)^(2/3) * exp(-2*x/3) * (x^2 - 3*x + 6) / (x - 3)^2:

`diff/mbrxc_x` := proc(Q, g)
  - (32*Pi)^(2/3)/6 * diff(Q, g)/(Q^2 * mbrxc_aux_dfdx(mbrxc_x(Q)))
end proc:

mbrxc_Q := (x, t) ->
        mbrxc_a1*(2*t) - K_FACTOR_C + mbrxc_a2*x^2 + mbrxc_a3*x^4:

(* This is unstable for small x *)
mbrxc_v0 := x ->
  - (32*Pi)^(1/3)/(8*X_FACTOR_C) * exp(x/3)*(8 - exp(-x)*(x^2 + 5*x + 8))/(x*(1 + x)^(1/3)):

(* so we use a series expansion for small x, remembering that we need
up to 4th derivatives which means that we need to increase the
truncation order by 4. *)
mbrxc_v_smallq := x -> eval(convert(taylor(mbrxc_v0(y), y = 0, 9), polynom), y=x):
mbrxc_v_cutoff := DBL_EPSILON^(1/4):
mbrxc_v := x -> my_piecewise3(x < mbrxc_v_cutoff, mbrxc_v_smallq(x), mbrxc_v0(m_max(x, mbrxc_v_cutoff))):

mbrxc_f := (x, u, t) ->
  - mbrxc_v(mbrxc_x(mbrxc_Q(x, t)))/2:

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(mbrxc_f, rs, z, xs0, xs1, u0, u1, t0, t1):
