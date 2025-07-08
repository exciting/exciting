(*
 Copyright (C) 2024 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

(* parameters from page 3 *)
lak_h0x := 1.174:
lak_mu_ax := -(97 + 3*lak_h0x + sqrt((3*lak_h0x)^2 + 74166*lak_h0x - 64175))/1200:
lak_nu_a := (73-50*lak_mu_ax)/5000:
lak_mu_sx := (10+60*lak_mu_ax)/81:
lak_nu_s := -(1606 - 50*lak_mu_ax)/18225:
lak_bx := 4.9479:
lak_ax := 1.1:
lak_anum := 5:

(* SI eq 2 *)
lak_fsa   := (s, a) -> lak_h0x*lak_gx(s) + (1-lak_fx(a))*(lak_h1x(s)-lak_h0x)*lak_gnum(s):
(* SI eq 4 *)
lak_gx := s -> 1 - m_recexp(sqrt(s)/lak_bx):

(* SI eq 5. This term has poor behavior around a=0, so we have to do a series expansion *)
lak_fx0 := a -> 2/Pi * arctan(Pi/2*(lak_c1*(a-1)/a + lak_c2*(a-1)^2)):
lak_fx_taylor := a -> eval(convert(eval(series(lak_fx0(b),b=0,9),csgn=-1),polynom),b=a):
lak_fx := a -> my_piecewise3(a <= DBL_EPSILON, lak_fx_taylor(a), lak_fx0(m_max(a,DBL_EPSILON))):
(* SI eq 6 *)
lak_c1 := lak_mu_ax/(lak_h0x-1):
(* SI eq 7 *)
lak_c2 := (lak_mu_ax + lak_nu_a)/(lak_h0x-1):
(* SI eq 8 *)
lak_h1x := s -> lak_hx_ge4(s) + lak_kx(s)*(lak_ax - lak_hx_ge4(s)):
(* SI eq 9 *)
lak_hx_ge4 := s -> 1 + lak_mu_sx*s^2 + lak_nu_s*s^4 + lak_h0x*(1-lak_gx(s)):
(* SI eq 10 *)
lak_kx := s -> m_recexp((s/lak_ax)^2 * (1+s^2)):
(* SI eq 11 *)
lak_gnum := s -> 1 - m_recexp((s/lak_anum)^2):

(* Build the functional *)
lak_alpha := (x, t) -> (t/K_FACTOR_C) * m_max(1 - x^2/(8*t), 1e-10):
lak_f   := (x, u, t) -> lak_fsa(x*X2S, lak_alpha(x,t)):
f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(lak_f, rs, z, xs0, xs1, u0, u1, t0, t1):
