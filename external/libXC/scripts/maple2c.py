#!/usr/bin/env python3

# Copyright (C) 2021 M.A.L. Marques
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import sys, re, os
from argparse import ArgumentParser
from maple2c_lib.utils import *

def maple2c_init():
  mparser = ArgumentParser(usage="Convert a maple file into C code")
  mparser.add_argument('--srcdir', type = str, default = ".",
                       help='Directory where to find the source code')
  mparser.add_argument('--functional', type = str,
                       help='Name of the functional', required = True)
  mparser.add_argument('--maxorder', type = int, default = 3,
                       help='Maximum order of the derivatives')
  mparser.add_argument('--simplify', type = str, default = False,
                       help='Should we add a final simplification command?')

  parse = mparser.parse_args()

  # gather all information necessary to maple2c
  params = {}
  params["srcdir"]     = parse.srcdir
  params["functional"] = parse.functional
  params["maxorder"]   = parse.maxorder
  params["simplify_begin"] = "simplify(" if parse.simplify else ""
  params["simplify_end"]   = ", symbolic)" if parse.simplify else ""

  # extract the family of the functional
  m = re.match(r'^(hyb_)?([^_]*)', params["functional"])
  params["family"] = m.group(2)

  # find out where the maple file resides
  possible_paths = [
    params["functional"],
    params["srcdir"] + "/maple/" + params["functional"],
    params["srcdir"] + "/maple/" + params["family"] + "_vxc/" + params["functional"],
    params["srcdir"] + "/maple/" + params["family"] + "_exc/" + params["functional"],
    ]
  for p in possible_paths:
    if os.path.isfile(p):
      params["maple_file"] = p
      break
    if os.path.isfile(p + ".mpl"):
      params["maple_file"] = p + ".mpl"
      break

  if "maple_file" not in params:
    print("File '" + params["functional"] + ".mpl" + "' not found")
    sys.exit(1)

  # we now read the header of the maple file to configure the functional
  fh = open(params["maple_file"], "r")
  params["replace"] = []
  params["prefix"]  = ""
  for line in fh:
    m = re.match(r'^\(\* type:\s(\S*)\s', line)
    if m:
      params["functype"] = m.group(1)

    m = re.match(r'^\(\* replace:\s*"([^"]*)"\s*->\s*"([^"]*)"', line)
    if m:
      params["replace"].append(m.group(1, 2))

    if re.match(r'^\(\* prefix:', line):
      for line in fh:
        if re.match(r'^\*', line):
          break
        else:
          params["prefix"] += line

  if "functype" not in params:
    print("Could not determine type of functional")
    print("Please add something like '(* type: lda_exc *)' to the maple file")
    sys.exit(1)

  return params

# the code starts here
params = maple2c_init()

from maple2c_lib.lda  import work_lda_exc, work_lda_vxc
from maple2c_lib.gga  import work_gga_exc, work_gga_vxc
from maple2c_lib.mgga import work_mgga_exc, work_mgga_vxc

if params["functype"] == "lda_exc":
  work_lda_exc(params)

elif params["functype"] == "lda_vxc":
  work_lda_vxc(params)

elif params["functype"] == "gga_exc":
  work_gga_exc(params)

elif params["functype"] == "gga_vxc":
  work_gga_vxc(params)

elif params["functype"] == "mgga_exc":
  work_mgga_exc(params)

elif params["functype"] == "mgga_vxc":
  work_mgga_vxc(params)
