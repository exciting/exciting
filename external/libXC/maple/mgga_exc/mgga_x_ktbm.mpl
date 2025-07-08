(* type: mgga_exc *)
(* prefix:
  mgga_x_ktbm_params *params;

  assert(p->params != NULL);
  params = (mgga_x_ktbm_params * )(p->params);
*)

ktbm_p := x -> X2S^2*x^2:
ktbm_t := t -> t/K_FACTOR_C:

ktbm_top := (p, t) -> params_a_ct + params_a_at*p + params_a_bt*t + params_a_a2t*p*p + params_a_b2t*t*t + params_a_xt*p*t :
ktbm_bot := (p, t) -> params_a_cb + params_a_ab*p + params_a_bb*t + params_a_a2b*p*p + params_a_b2b*t*t + params_a_xb*p*t :

ktbm_f := (x, u, t) -> ktbm_top(ktbm_p(x),ktbm_t(t))/ktbm_bot(ktbm_p(x),ktbm_t(t)) :

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(ktbm_f, rs, z, xs0, xs1, u0, u1, t0, t1):
