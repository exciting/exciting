#!/usr/bin/env python3

# Copyright (C) 2021 M.A.L. Marques
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from maple2c_lib.utils import *

# these are the variables that the functional depends on
variables = ["rho_0_", "rho_1_", "sigma_0_", "sigma_1_", "sigma_2_"]

# get arguments of the functions
input_args  = "const double *rho, const double *sigma"
output_args = "xc_gga_out_params *out"

# the definition of the derivatives that libxc transmits to the calling program
partials = [
  ["zk"],
  ["vrho", "vsigma"],
  ["v2rho2", "v2rhosigma", "v2sigma2"],
  ["v3rho3", "v3rho2sigma", "v3rhosigma2", "v3sigma3"],
  ["v4rho4", "v4rho3sigma", "v4rho2sigma2", "v4rhosigma3", "v4sigma4"]
]

#####################################################################
def work_gga_exc(params):
  '''Process a GGA functional for the energy'''

  derivatives = partials_to_derivatives(params, "gga", partials)
  
  der_def, out_c = maple_define_derivatives(variables, derivatives, "mf")
  
  out_c = ", ".join(out_c)
  if out_c != "": out_c = ", " + out_c

  # we join all the pieces
  maple_code  = '''
# zk is energy per unit particle
mzk  := (r0, r1, s0, s1, s2) -> \\
  {} + \\
    f(r_ws(dens(r0, r1)), zeta(r0, r1), xt(r0, r1, s0, s1, s2), xs0(r0, r1, s0, s2), xs1(r0, r1, s0, s2)) \\
  {} :

  (* mf is energy per unit volume *)
  mf   := (r0, r1, s0, s1, s2) -> eval(dens(r0, r1)*mzk(r0, r1, s0, s1, s2)):

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

{}

{}
C([{}{}], optimize, deducetypes=false):

'''.format(der_def, maple_code, maple_zk, out_c)
  }

  maple2c_run(params, variables, derivatives, variants, 0, input_args, output_args)


#####################################################################
def work_gga_vxc(params):
  '''Process a GGA functional for the potential'''

  all_derivatives = partials_to_derivatives(params, "gga", partials)

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
mzk  := (r0, r1, s0, s1, s2) -> \\
  {} + \\
    f(r_ws(dens(r0, r1)), zeta(r0, r1), xt(r0, r1, s0, s1, s2), xs0(r0, r1, s0, s2), xs1(r0, r1, s0, s2)) \\
  {} :

(* mf is the up potential *)
mf0   := (r0, r1, s0, s1, s2) -> eval(mzk(r0, r1, s0, s1, s2)):
mf1   := (r0, r1, s0, s1, s2) -> eval(mzk(r1, r0, s2, s1, s0)):

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

{}

{}
C([{}, {}{}], optimize, deducetypes=false):

'''.format(der_def_pol, maple_code, maple_vrho0, maple_vrho1, out_c_pol)
  }

  maple2c_run(params, variables, derivatives, variants, 1, input_args, output_args)
