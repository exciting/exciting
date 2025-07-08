#!/usr/bin/env python3

# Copyright (C) 2021 M.A.L. Marques
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from maple2c_lib.utils import *

# these are the variables that the functional depends on
variables = ["rho_0_", "rho_1_", "sigma_0_", "sigma_1_", "sigma_2_", "lapl_0_", "lapl_1_", "tau_0_", "tau_1_"]
  
# get arguments of the functions
input_args  = "const double *rho, const double *sigma, const double *lapl, const double *tau"
output_args = "xc_mgga_out_params *out"

# the definition of the derivatives that libxc transmits to the calling program
partials = [
  ["zk"],
  ["vrho", "vsigma", "vlapl", "vtau"],
  ["v2rho2", "v2rhosigma", "v2rholapl", "v2rhotau", "v2sigma2",
   "v2sigmalapl", "v2sigmatau", "v2lapl2", "v2lapltau", "v2tau2"],
  ["v3rho3", "v3rho2sigma", "v3rho2lapl", "v3rho2tau", "v3rhosigma2",
   "v3rhosigmalapl", "v3rhosigmatau", "v3rholapl2", "v3rholapltau",
   "v3rhotau2", "v3sigma3", "v3sigma2lapl", "v3sigma2tau", "v3sigmalapl2",
   "v3sigmalapltau", "v3sigmatau2", "v3lapl3", "v3lapl2tau", "v3lapltau2",
   "v3tau3"],
  ["v4rho4", "v4rho3sigma", "v4rho3lapl", "v4rho3tau", "v4rho2sigma2",
   "v4rho2sigmalapl", "v4rho2sigmatau", "v4rho2lapl2", "v4rho2lapltau",
   "v4rho2tau2", "v4rhosigma3", "v4rhosigma2lapl", "v4rhosigma2tau",
   "v4rhosigmalapl2", "v4rhosigmalapltau", "v4rhosigmatau2",
   "v4rholapl3", "v4rholapl2tau", "v4rholapltau2", "v4rhotau3",
   "v4sigma4", "v4sigma3lapl", "v4sigma3tau", "v4sigma2lapl2",
   "v4sigma2lapltau", "v4sigma2tau2", "v4sigmalapl3", "v4sigmalapl2tau",
   "v4sigmalapltau2", "v4sigmatau3", "v4lapl4", "v4lapl3tau",
   "v4lapl2tau2", "v4lapltau3", "v4tau4",
  ]
]

#####################################################################
def work_mgga_exc(params):
  '''Process a MGGA functional for the energy'''

  derivatives = partials_to_derivatives(params, "mgga", partials)
  
  der_def, out_c = maple_define_derivatives(variables, derivatives, "mf")
  
  out_c = ", ".join(out_c)
  if out_c != "": out_c = ", " + out_c

  # we join all the pieces
  maple_code  = '''
# zk is energy per unit particle
mzk  := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> \\
  {} + \\
    f(r_ws(dens(r0, r1)), zeta(r0, r1), xt(r0, r1, s0, s1, s2), xs0(r0, r1, s0, s2), xs1(r0, r1, s0, s2), u0(r0, r1, l0, l1), u1(r0, r1, l0, l1), t0(r0, r1, tau0, tau1), t1(r0, r1, tau0, tau1)) \\
  {} :

  (* mf is energy per unit volume *)
  mf   := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> eval(dens(r0, r1)*mzk(r0, r1, s0, s1, s2, l0, l1, tau0, tau1)):

$include <util.mpl>
'''.format(params["simplify_begin"], params["simplify_end"])
  
  maple_zk = " zk_0_ = mzk(" + ", ".join(variables) + ")"

  # we build 2 variants of the functional, for unpolarized, and polarized densities
  variants = {
    "unpol": '''
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):

{}

{}
C([{}{}], optimize, deducetypes=false):

'''.format(der_def, maple_code, maple_zk, out_c),

    "pol": '''
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma2)/r1^(1 + 1/DIMENSIONS):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0 + 2*sigma1 + sigma2)/(r0 + r1)^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> l0/(r0^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> l1/(r1^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> tau0/(r0^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> tau1/(r1^(1 + 2/DIMENSIONS)):

{}

{}
C([{}{}], optimize, deducetypes=false):

'''.format(der_def, maple_code, maple_zk, out_c)
  }

  maple2c_run(params, variables, derivatives, variants, 0, input_args, output_args)


#####################################################################
def work_mgga_vxc(params):
  '''Process a MGGA functional for the potential'''

  all_derivatives = partials_to_derivatives(params, "mgga", partials)

  derivatives, derivatives1, derivatives2 = filter_vxc_derivatives(all_derivatives)
  
  # we obtain the missing pieces for maple
  # unpolarized calculation
  der_def_unpol, out_c_unpol = maple_define_derivatives(variables, derivatives1, "mf0")
  out_c_unpol = ", ".join(out_c_unpol)
  if out_c_unpol != "": out_c_unpol = ", " + out_c_unpol

  # polarized calculation
  der_def_pol1, out_c_pol1 = maple_define_derivatives(variables, derivatives1, "mf0")
  der_def_pol2, out_c_pol2 = maple_define_derivatives(variables, derivatives2, "mf1")
  
  der_def_pol = der_def_pol1 + der_def_pol2
  out_c_pol   = ", ".join(sorted(out_c_pol1 + out_c_pol2, key=sort_alphanumerically))
  if out_c_pol != "": out_c_pol = ", " + out_c_pol
  
  # we join all the pieces
  maple_code  = '''
mzk  := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> \\
  {} + \\
    f(r_ws(dens(r0, r1)), zeta(r0, r1), xt(r0, r1, s0, s1, s2), xs0(r0, r1, s0, s2), xs1(r0, r1, s0, s2), u0(r0, r1, l0, l1), u1(r0, r1, l0, l1), t0(r0, r1, tau0, tau1), t1(r0, r1, tau0, tau1)) \\
  {} :

(* mf is the up potential *)
mf0   := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> eval(mzk(r0, r1, s0, s1, s2, l0, l1, tau0, tau1)):
mf1   := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> eval(mzk(r1, r0, s2, s1, s0, l1, l0, tau1, tau0)):

$include <util.mpl>
'''.format(params["simplify_begin"], params["simplify_end"])
  
  maple_vrho0 = " vrho_0_ = mf0(" + ", ".join(variables) + ")"
  maple_vrho1 = " vrho_1_ = mf1(" + ", ".join(variables) + ")"

  # we build 2 variants of the functional, for unpolarized, and polarized densities
  variants = {
    "unpol": '''
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):

{}

{}
C([{}{}], optimize, deducetypes=false):

'''.format(der_def_unpol, maple_code, maple_vrho0, out_c_unpol),

    "pol": '''
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma2)/r1^(1 + 1/DIMENSIONS):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0 + 2*sigma1 + sigma2)/(r0 + r1)^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> l0/(r0^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> l1/(r1^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> tau0/(r0^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> tau1/(r1^(1 + 2/DIMENSIONS)):

{}

{}
C([{}, {}{}], optimize, deducetypes=false):

'''.format(der_def_pol, maple_code, maple_vrho0, maple_vrho1, out_c_pol)
  }

  maple2c_run(params, variables, derivatives, variants, 1, input_args, output_args)
