(*
 2018 Authored by Andrea Kreppel
 2022 Edited by Henryk Laqua

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 Short-range PBE exchange functional Goll/Werner/Stoll
 E. Goll, H.-J. Werner, and H. Stoll., Phys. Chem. Chem. Phys. 7, 3917 (2005).
 DOI:10.1039/B509242F
*)

(* type: gga_exc *)

(* short-range LDA is the basis of short-range PBE*)
$include "lda_x_erf.mpl"

(*
 in short-range PBE, the constant b is now mu-dependent.
*)

(*coefficients given in the text below eq. (5)*)
pbe_x_erf_gws_c1_b := x -> 1 + 22*x^2 + 144*x^4:
pbe_x_erf_gws_c2_b := x -> 2*x^2*(-7+72*x^2):
pbe_x_erf_gws_c3_b := x -> -864*x^4*(-1+2*x^2):
pbe_x_erf_gws_c4_b := x -> x^2*(-3 - 24*x^2 + 32*x^4 + 8*x*sqrt(Pi)*erf(1/(2*x))):
(*exponential term of eq. (4)*)
exp_b := mu_t -> exp(1/(4*(mu_t)^2)):

(*eq. (4)*)
pbe_x_erf_gws_b := mu_t -> ( -pbe_x_erf_gws_c1_b(mu_t) + pbe_x_erf_gws_c2_b(mu_t)*exp_b(mu_t) )/( pbe_x_erf_gws_c3_b(mu_t) + 54*pbe_x_erf_gws_c4_b(mu_t)*exp_b(mu_t) ):

#(*Taylor expansion around 0 to avoid numerical problems with small mu_t*)
#pbe_x_erf_gws_b_small_mu := mu_t -> 0.8641975309e-1+.1210182562*mu_t+0.6151325935e-1*mu_t^2-.2514737006*mu_t^3-1.756502351*mu_t^4-7.573244644*mu_t^5:
(*the limit mu_t -> 0 is dominated by exp_b*)
pbe_x_erf_gws_b_small_mu := mu_t -> pbe_x_erf_gws_c2_b(mu_t) /(54*pbe_x_erf_gws_c4_b(mu_t)):

(*Expansion around infinity to avoid numerical problems with large mu_t*)
pbe_x_erf_gws_b_large_mu := mu_t -> 1/(72*mu_t^2) - 1/(17280*mu_t^4) - 23/(358400*mu_t^6):

(*
 switch between Taylor expansions and exact form with differentiable step function
*)
# exp_b explodes otherwise (exp_b(0.5) = 2.7e+43) (we also want some spare-room for derivatives)
pbe_x_erf_gws_b_thresh_small := 0.05:
# Should be fine. largest scaling is O(mu_t^4)=O(1e40) (we also want some spare-room for derivatives)
pbe_x_erf_gws_b_thresh_large := 1e10:
pbe_x_erf_gws_b_piece := mu_t -> my_piecewise5(mu_t < pbe_x_erf_gws_b_thresh_small, pbe_x_erf_gws_b_small_mu(mu_t),
                                          mu_t > pbe_x_erf_gws_b_thresh_large, pbe_x_erf_gws_b_large_mu(mu_t),
                                          pbe_x_erf_gws_b(mu_t)):

(*
 special case b(0) == 7/81
*)
pbe_x_erf_gws_b_piece0 := 7/81:


(*fixed (not mu, rs,z dependent) PBE parameters*)
(* prefix:
  gga_x_pbe_erf_gws_params *params;
  assert(p->params != NULL);
  params = (gga_x_pbe_erf_gws_params * )(p->params);
*)

(*default parameters*)
$ifdef gga_x_pbe_erf_gws_params
params_a_kappa := 0.8040:
params_a_b_PBE := 0.2195149727645171:
params_a_ax    := 19.0:
$endif

pbe_x_erf_gws_kappa_fx := (rs,z) -> params_a_kappa: (* we may want to extent this to a version with a density-depended kappa later *)
pbe_x_erf_gws_x_b_orig := params_a_b_PBE:
pbe_x_erf_gws_ax := params_a_ax:

(*modified b from eq. (5)*)
pbe_x_erf_gws_b_mod := mu_t -> pbe_x_erf_gws_x_b_orig/pbe_x_erf_gws_b_piece0 * pbe_x_erf_gws_b_piece(mu_t) * exp(-pbe_x_erf_gws_ax*mu_t^2):

#we emplow nu/2 (spin-nu) instead of nu 
nu_2 := (rs,z) -> nu(rs,z)/2:

(* second part of eq. (3)*)
pbe_x_erf_gws_Fx := (rs,z,s) -> 1 + pbe_x_erf_gws_kappa_fx(rs,z)*(1 - pbe_x_erf_gws_kappa_fx(rs,z)/(pbe_x_erf_gws_kappa_fx(rs,z) + pbe_x_erf_gws_b_mod(nu_2(rs,z))*s^2)):

$include "lda_x_erf.mpl"
(* first part of eq. (3)*)
f_pbe_x_erf_gws_spin := (rs, z, xs) -> lda_x_erf_spin(rs,z)*pbe_x_erf_gws_Fx(rs,z,xs*X2S):

rs_a := (rs,z) -> simplify(r_ws(n_spin(rs,z))):
rs_b := (rs,z) -> simplify(r_ws(n_spin(rs,-z))):
f := (rs, z, xt, xs0, xs1) -> simplify((
+ my_piecewise3(screen_dens(rs, z),0,f_pbe_x_erf_gws_spin(rs_a(rs,z),1,xs0)*n_spin(rs, z))
+ my_piecewise3(screen_dens(rs,-z),0,f_pbe_x_erf_gws_spin(rs_b(rs,z),1,xs1)*n_spin(rs,-z))
)/n_total(rs)):
