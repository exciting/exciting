#!/usr/bin/env python3

# Copyright (C) 2021 M.A.L. Marques
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import sys, os, re
import subprocess
from textwrap import wrap
from maple2c_lib.utils import enumerate_spin_partials

replace = [
  (r"\[n0, s0, u0\]", ""),
  (r"\[n1, s2, u1\]", ""),
  (r"\[n0, n1, s0, s1, s2, u0, u1, ked1, ked2\]" , ""),
  (r"Derivative", ""),

  (r"\[1, 0, 0\]\[ked1\]", "ked1_vrho[0]"),
  (r"\[0, 1, 0\]\[ked1\]", "ked1_vsigma[0]"),
  (r"\[0, 0, 1\]\[ked1\]", "ked1_vlapl[0]"),

  (r"\[2, 0, 0\]\[ked1\]", "ked1_v2rho2[0]"),
  (r"\[1, 1, 0\]\[ked1\]", "ked1_v2rhosigma[0]"),
  (r"\[1, 0, 1\]\[ked1\]", "ked1_v2rholapl[0]"),
  (r"\[0, 2, 0\]\[ked1\]", "ked1_v2sigma2[0]"),
  (r"\[0, 1, 1\]\[ked1\]", "ked1_v2sigmalapl[0]"),
  (r"\[0, 0, 2\]\[ked1\]", "ked1_v2lapl2[0]"),

  (r"\[3, 0, 0\]\[ked1\]", "ked1_v3rho3[0]"),
  (r"\[2, 1, 0\]\[ked1\]", "ked1_v3rho2sigma[0]"),
  (r"\[2, 0, 1\]\[ked1\]", "ked1_v3rho2lapl[0]"),
  (r"\[1, 2, 0\]\[ked1\]", "ked1_v3rhosigma2[0]"),
  (r"\[1, 1, 1\]\[ked1\]", "ked1_v3rhosigmalapl[0]"),
  (r"\[1, 0, 2\]\[ked1\]", "ked1_v3rholapl2[0]"),
  (r"\[0, 3, 0\]\[ked1\]", "ked1_v3sigma3[0]"),
  (r"\[0, 2, 1\]\[ked1\]", "ked1_v3sigma2lapl[0]"),
  (r"\[0, 1, 2\]\[ked1\]", "ked1_v3sigmalapl2[0]"),
  (r"\[0, 0, 3\]\[ked1\]", "ked1_v3lapl3[0]"),

  (r"\[4, 0, 0\]\[ked1\]", "ked1_v4rho4[0]"),
  (r"\[3, 1, 0\]\[ked1\]", "ked1_v4rho3sigma[0]"),
  (r"\[3, 0, 1\]\[ked1\]", "ked1_v4rho3lapl[0]"),
  (r"\[2, 2, 0\]\[ked1\]", "ked1_v4rho2sigma2[0]"),
  (r"\[2, 1, 1\]\[ked1\]", "ked1_v4rho2sigmalapl[0]"),
  (r"\[2, 0, 2\]\[ked1\]", "ked1_v4rho2lapl2[0]"),
  (r"\[1, 3, 0\]\[ked1\]", "ked1_v4rhosigma3[0]"),
  (r"\[1, 2, 1\]\[ked1\]", "ked1_v4rhosigma2lapl[0]"),
  (r"\[1, 1, 2\]\[ked1\]", "ked1_v4rhosigmalapl2[0]"),
  (r"\[1, 0, 3\]\[ked1\]", "ked1_v4rholapl3[0]"),
  (r"\[0, 4, 0\]\[ked1\]", "ked1_v4sigma4[0]"),
  (r"\[0, 3, 1\]\[ked1\]", "ked1_v4sigma3lapl[0]"),
  (r"\[0, 2, 2\]\[ked1\]", "ked1_v4sigma2lapl2[0]"),
  (r"\[0, 1, 3\]\[ked1\]", "ked1_v4sigmalapl3[0]"),
  (r"\[0, 0, 4\]\[ked1\]", "ked1_v4lapl4[0]"),
]

ders = (
  "vrho", "vsigma", "vlapl", "vtau",

  "v2rho2", "v2rhosigma", "v2rholapl", "v2rhotau", "v2sigma2",
  "v2sigmalapl", "v2sigmatau", "v2lapl2", "v2lapltau", "v2tau2",

  "v3rho3", "v3rho2sigma", "v3rho2lapl", "v3rho2tau",
  "v3rhosigma2", "v3rhosigmalapl", "v3rhosigmatau", "v3rholapl2",
  "v3rholapltau", "v3rhotau2", "v3sigma3", "v3sigma2lapl",
  "v3sigma2tau", "v3sigmalapl2", "v3sigmalapltau", "v3sigmatau2",
  "v3lapl3", "v3lapl2tau", "v3lapltau2", "v3tau3",

  "v4rho4", "v4rho3sigma", "v4rho3lapl", "v4rho3tau",
  "v4rho2sigma2", "v4rho2sigmalapl", "v4rho2sigmatau", "v4rho2lapl2",
  "v4rho2lapltau", "v4rho2tau2", "v4rhosigma3", "v4rhosigma2lapl",
  "v4rhosigma2tau", "v4rhosigmalapl2", "v4rhosigmalapltau", "v4rhosigmatau2",
  "v4rholapl3", "v4rholapl2tau", "v4rholapltau2", "v4rhotau3",
  "v4sigma4", "v4sigma3lapl", "v4sigma3tau", "v4sigma2lapl2",
  "v4sigma2lapltau", "v4sigma2tau2", "v4sigmalapl3", "v4sigmalapl2tau",
  "v4sigmalapltau2", "v4sigmatau3", "v4lapl4", "v4lapl3tau",
  "v4lapl2tau2", "v4lapltau3", "v4tau4"
)

params = {}
params["maxorder"] = 4

ders_def = []
for der in ders:
  ders_def.extend(enumerate_spin_partials(der, "mgga"))

# create replaces. E.g.
#  ('\[1, 0, 0, 0, 0, 0, 0, 0, 0\]\[mgga\]' , "mgga_vrho[0]")
for der in ders_def:
  der[1] = re.sub(r"_(\d+)_", r"[\1]", der[1])
  replace.append(
    (r"\[" + ", ".join([str(i) for i in der[0]]) + r"\]\[mgga\]",
     "mgga_" + str(der[1])))

out1 = ""
out2 = ""
for der in ders_def:
  # let us do one order each time
  if not re.match(r"v" + sys.argv[1] + "[a-z]", der[1]): # argv[1] = "", "2", "3", or "4"
    continue

  # no need to calculate tau derivatives
  if re.search(r"tau", der[1]): # v2, v3, v4
    continue
  
  mvars = ("n0", "n1", "s0", "s1", "s2", "u0", "u1", "t0", "t1")

  mder = ",".join(["{{{}, {}}}".format(mvars[i], der[0][i]) for i in range(len(mvars))])

  # run mathematica
  fh = open("/tmp/math.m", "w")
  fh.write("Print[ToString[FullSimplify[D[mgga[n0, n1, s0, s1, s2, u0, u1, ked1[n0, s0, u0], ked2[n1, s2, u1]], " + mder + "]], FormatType -> InputForm]]")
  fh.close()

  run = subprocess.run(
    ["math", "-script", "/tmp/math.m"],
    stdout=subprocess.PIPE, universal_newlines=True)
  out = run.stdout.strip()
  
  # make the appropriate replacements in mathematica output
  for str1, str2 in replace:
    out = re.sub(str1, str2, out)

    # repeat the substitutions for ked2
    str1 = re.sub(r"ked1", "ked2", str1)
    str2 = re.sub(r"ked1", "ked2", str2)
    out = re.sub(str1, str2, out)

  # let us get rid of the powers
  out = re.sub(r"([A-Za-z][A-Za-z_\[\]0-9]*)\^2", r"\1*\1",     out)
  out = re.sub(r"([A-Za-z][A-Za-z_\[\]0-9]*)\^3", r"\1*\1*\1",   out)
  out = re.sub(r"([A-Za-z][A-Za-z_\[\]0-9]*)\^4", r"\1*\1*\1*\1", out)

  out = der[1] + " = " + out + ";\n"

  res = re.search(r"(.\d?)(.*?)\[(\d+)\]", der[1])

  if res and res.group(3) == "0":
    out1 += "\n".join(wrap(out, initial_indent='', subsequent_indent='  ',
                           expand_tabs=True), width=100) + "\n"
  else:  
    out2 += "\n".join(wrap(out, initial_indent='  ', subsequent_indent='    ',
                           expand_tabs=True), width=100) + "\n"

print(out1 + "\nif(func->nspin == XC_POLARIZED){\n" + out2 + "}\n")
