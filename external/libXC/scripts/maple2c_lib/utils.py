#!/usr/bin/env python3

# Copyright (C) 2021 M.A.L. Marques
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import sys, os, re, subprocess

# we need this in a couple of places
der_name = ("EXC", "VXC", "FXC", "KXC", "LXC", "MXC")

#####################################################################
# sort by character and then by number
def sort_alphanumerically(x):
  res = re.match(r"([^_]+)_([0-9]+)_", x)
  return (res.group(1), int(res.group(2)))


#####################################################################
def partials_to_derivatives(params, func_type, partials):
  derivatives = []
  for order in range(params["maxorder"] + 1):
    derivatives.append([])
    for der in partials[order]:
      derivatives[order].extend(enumerate_spin_partials(der, func_type))

  return derivatives


#####################################################################
def filter_vxc_derivatives(all_derivatives):
  '''This separates the derivatives of vxc into derivatives of vxc_0
  and vxc_1. All other derivatives (e.g. vsigma) are ignored.'''

  derivatives  = []
  derivatives1 = []
  derivatives2 = []

  for order in range(len(all_derivatives) - 1):
    derivatives.append([])
    derivatives1.append([])
    derivatives2.append([])
    
    for der in all_derivatives[order + 1]:
      if der[0][0] > 0:
        der[0][0] -= 1
        derivatives1[order].append(der)
        derivatives [order].append(der)
      elif der[0][1] > 0:
        der[0][1] -= 1
        derivatives2[order].append(der)
        derivatives [order].append(der)

  return derivatives, derivatives1, derivatives2


def enumerate_spin_partials(derivative, func_type):
  '''Given the name of a derivative (such as 'v2rho2')
  and a functional type ("lda", "gga", "mgga"), return all spin 
  variants such as

  [
    [[2, 0], 'v2rho2_0_'], 
    [[1, 1], 'v2rho2_1_'], 
    [[0, 2], 'v2rho2_2_']
  ]
  '''

  words = ("rho", "sigma", "lapl", "tau")
  partials = {
    "rho"   :  [[[0, 0, 0, 0, 0, 0, 0, 0, 0]],    # 0th-order
                [[1, 0, 0, 0, 0, 0, 0, 0, 0],     # 1st-order
                 [0, 1, 0, 0, 0, 0, 0, 0, 0]],
                [[2, 0, 0, 0, 0, 0, 0, 0, 0],     # 2nd-order
                 [1, 1, 0, 0, 0, 0, 0, 0, 0],
                 [0, 2, 0, 0, 0, 0, 0, 0, 0]],
                [[3, 0, 0, 0, 0, 0, 0, 0, 0],     # 3rd-order
                 [2, 1, 0, 0, 0, 0, 0, 0, 0],
                 [1, 2, 0, 0, 0, 0, 0, 0, 0],
                 [0, 3, 0, 0, 0, 0, 0, 0, 0]],
                [[4, 0, 0, 0, 0, 0, 0, 0, 0],     # 4th-order
                 [3, 1, 0, 0, 0, 0, 0, 0, 0],
                 [2, 2, 0, 0, 0, 0, 0, 0, 0],
                 [1, 3, 0, 0, 0, 0, 0, 0, 0],
                 [0, 4, 0, 0, 0, 0, 0, 0, 0]],
    ],
    "sigma" :  [[[0, 0, 0, 0, 0, 0, 0, 0, 0]],    # 0th-order
                [[0, 0, 1, 0, 0, 0, 0, 0, 0],     # 1st-order
                 [0, 0, 0, 1, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 1, 0, 0, 0, 0]],
                [[0, 0, 2, 0, 0, 0, 0, 0, 0],     # 2nd-order
                 [0, 0, 1, 1, 0, 0, 0, 0, 0],
                 [0, 0, 1, 0, 1, 0, 0, 0, 0],
                 [0, 0, 0, 2, 0, 0, 0, 0, 0],
                 [0, 0, 0, 1, 1, 0, 0, 0, 0],
                 [0, 0, 0, 0, 2, 0, 0, 0, 0]],
                [[0, 0, 3, 0, 0, 0, 0, 0, 0],     # 3rd-order
                 [0, 0, 2, 1, 0, 0, 0, 0, 0],
                 [0, 0, 2, 0, 1, 0, 0, 0, 0],
                 [0, 0, 1, 2, 0, 0, 0, 0, 0],
                 [0, 0, 1, 1, 1, 0, 0, 0, 0],
                 [0, 0, 1, 0, 2, 0, 0, 0, 0],
                 [0, 0, 0, 3, 0, 0, 0, 0, 0],
                 [0, 0, 0, 2, 1, 0, 0, 0, 0],
                 [0, 0, 0, 1, 2, 0, 0, 0, 0],
                 [0, 0, 0, 0, 3, 0, 0, 0, 0]],
                [[0, 0, 4, 0, 0, 0, 0, 0, 0],     # 4th-order
                 [0, 0, 3, 1, 0, 0, 0, 0, 0],
                 [0, 0, 3, 0, 1, 0, 0, 0, 0],
                 [0, 0, 2, 2, 0, 0, 0, 0, 0],
                 [0, 0, 2, 1, 1, 0, 0, 0, 0],
                 [0, 0, 2, 0, 2, 0, 0, 0, 0],
                 [0, 0, 1, 3, 0, 0, 0, 0, 0],
                 [0, 0, 1, 2, 1, 0, 0, 0, 0],
                 [0, 0, 1, 1, 2, 0, 0, 0, 0],
                 [0, 0, 1, 0, 3, 0, 0, 0, 0],
                 [0, 0, 0, 4, 0, 0, 0, 0, 0],
                 [0, 0, 0, 3, 1, 0, 0, 0, 0],
                 [0, 0, 0, 2, 2, 0, 0, 0, 0],
                 [0, 0, 0, 1, 3, 0, 0, 0, 0],
                 [0, 0, 0, 0, 4, 0, 0, 0, 0]]
    ],
    "lapl" :   [[[0, 0, 0, 0, 0, 0, 0, 0, 0]],    # 0th-order
                [[0, 0, 0, 0, 0, 1, 0, 0, 0],     # 1st-order
                 [0, 0, 0, 0, 0, 0, 1, 0, 0]],
                [[0, 0, 0, 0, 0, 2, 0, 0, 0],     # 2nd-order
                 [0, 0, 0, 0, 0, 1, 1, 0, 0],
                 [0, 0, 0, 0, 0, 0, 2, 0, 0]],
                [[0, 0, 0, 0, 0, 3, 0, 0, 0],     # 3rd-order
                 [0, 0, 0, 0, 0, 2, 1, 0, 0],
                 [0, 0, 0, 0, 0, 1, 2, 0, 0],
                 [0, 0, 0, 0, 0, 0, 3, 0, 0]],
                [[0, 0, 0, 0, 0, 4, 0, 0, 0],     # 4th-order
                 [0, 0, 0, 0, 0, 3, 1, 0, 0],
                 [0, 0, 0, 0, 0, 2, 2, 0, 0],
                 [0, 0, 0, 0, 0, 1, 3, 0, 0],
                 [0, 0, 0, 0, 0, 0, 4, 0, 0]]
    ],
    "tau"  :   [[[0, 0, 0, 0, 0, 0, 0, 0, 0]],    # 0th-order
                [[0, 0, 0, 0, 0, 0, 0, 1, 0],     # 1st-order
                 [0, 0, 0, 0, 0, 0, 0, 0, 1]],
                [[0, 0, 0, 0, 0, 0, 0, 2, 0],     # 2nd-order
                 [0, 0, 0, 0, 0, 0, 0, 1, 1],
                 [0, 0, 0, 0, 0, 0, 0, 0, 2]],
                [[0, 0, 0, 0, 0, 0, 0, 3, 0],     # 3rd-order
                 [0, 0, 0, 0, 0, 0, 0, 2, 1],
                 [0, 0, 0, 0, 0, 0, 0, 1, 2],
                 [0, 0, 0, 0, 0, 0, 0, 0, 3]],
                [[0, 0, 0, 0, 0, 0, 0, 4, 0],     # 4th-order
                 [0, 0, 0, 0, 0, 0, 0, 3, 1],
                 [0, 0, 0, 0, 0, 0, 0, 2, 2],
                 [0, 0, 0, 0, 0, 0, 0, 1, 3],
                 [0, 0, 0, 0, 0, 0, 0, 0, 4]]
                ]
  }

  # finds out the order of each partial
  order = {}
  for word in words:
    order[word] = 0
  
    m = re.match(r'.*' + word + r'([0-9]*)', derivative)
    if m is not None:
      order[word] = 1 if m.group(1) == "" else int(m.group(1))

  # and this is the order of the derivative
  total_order = sum(order.values())

  max_n = {"lda": 2, "gga": 5, "mgga": 9}[func_type]

  all_derivatives = []
  der_n = 0
  for n_rho, p_rho in enumerate(partials["rho"][order["rho"]]):
    for n_sigma, p_sigma in enumerate(partials["sigma"][order["sigma"]]):
      for n_lapl, p_lapl in enumerate(partials["lapl"][order["lapl"]]):
        for n_tau, p_tau in enumerate(partials["tau"][order["tau"]]):
          # sum orders in all variables
            
          final_der = [0] * max_n
          for i in range(max_n):
            final_der[i] += \
                partials["rho"][order["rho"]][n_rho][i] + \
                partials["sigma"][order["sigma"]][n_sigma][i] + \
                partials["lapl"][order["lapl"]][n_lapl][i] + \
                partials["tau"][order["tau"]][n_tau][i]

          all_derivatives.append([final_der, derivative + "_" + str(der_n) + "_"])
          der_n += 1

  return all_derivatives


def maple_define_derivatives(variables, derivatives, func):
  '''Generates maple code to define derivatives. Output is

  out_derivatives = "
dmfd10 := (v0, v1) ->  eval(diff(mf(v0, v1), v0)):

dmfd01 := (v0, v1) ->  eval(diff(mf(v0, v1), v1)):
...
"
  out_cgeneration = "[
    'vrho_0_ = dmfd10(rho_0_, rho_1_)', 
    'vrho_1_ = dmfd01(rho_0_, rho_1_)'
    ...
  ]
  
  '''
  
  out_derivatives = ""
  out_cgeneration = []

  realvars = ", ".join(variables)

  for der_order in derivatives:
    for der in der_order:
      order, name = der
      to_derive = order.copy()
      
      # is there something to do?
      if all([v == 0 for v in to_derive]):
        break

      # derivates are always defined as a first derivative of another derivative
      # so we have to find out what is this previous function to derive
      for i_to_derive in range(len(to_derive)):
        if to_derive[i_to_derive] != 0:
          to_derive[i_to_derive] -= 1
          break

      if all([v == 0 for v in to_derive]):
        f_to_derive = func
      else:
        f_to_derive = "d" + func + "d" + "".join(str(i) for i in to_derive)

      # we build the expression of the derivative
      varname  = "d" + func + "d" + "".join(str(i) for i in order)
      vars     = "v" + ", v".join(str(i) for i in range(len(order)))
      derorder = ", ".join(str(i) for i in order)

      out_derivatives += varname + " := (" + vars + ") ->  eval(diff(" + \
        f_to_derive + "(" + vars +"), v" + str(i_to_derive)  + ")):\n\n"

      out_cgeneration.append(name + " = " + varname + "(" + realvars + ")")

  return out_derivatives, out_cgeneration


def print_c_header(params, out):
  cmd = "echo -e 'quit;' | maple 2>&1 | head -n 1 | sed 's/^.*Maple/Maple/'"
  maple_version = subprocess.check_output(cmd, shell=True).strip()

  out.write('''/*
  This file was generated automatically with {}.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : {}
  Maple source      : {}
  Type of functional: {}
*/

#define maple2c_order {}
'''.format(sys.argv[0], maple_version.decode(), params['maple_file'], params['functype'], params['maxorder']))


def maple2c_replace(text, extra_replace=()):
  '''Performs a series of string replacements in the maple generated C code'''
  
  # The replacements have to be made in order
  math_replace = (
    (r"_s_",     r"*"),
    (r"_a_",     r"->"),
    (r"_d_",     r"."),
    (r"_(\d+)_",  r"[\1]"),
    # convert constants like 0.225000000e-1 to 0.225e-1
    (r"0+e", r"e"),
    # convert numerical value of pi to constant
    (r"0.31415926535897932385e1", r"M_PI"),
    # have to do it here, as both Dirac(x) and Dirac(n, x) can appear
    (r"Dirac\(.*?\)", r"0.0"),
    # the derivative of the signum is 0 for us
    (r"signum\(1.*\)", r"0.0"),
    # optimizing specific calls to pow (sometimes we need to be careful with nested parenthesis)
    (r"pow\(0.1e1, (?:\(.*?\)|[^\(])*?\)",        r"0.1e1"),
    (r"pow\((.*?), *0.10*e1\)",                   r"(\1)"),
    (r"pow\((.*?), *-0.10*e1\)",                  r"0.1e1/(\1)"),
    (r"pow\((.*?), *0.15e1\)",                    r"POW_3_2(\1)"),
    (r"pow\((.*?), *-0.15e1\)",                   r"0.1e1/POW_3_2(\1)"),
    (r"pow\((.*?), *0.20*e1\)",                   r"POW_2(\1)"),
    (r"pow\((.*?), *-0.20*e1\)",                  r"0.1e1/POW_2(\1)"),
    (r"pow\((.*?), *0.30*e1\)",                   r"POW_3(\1)"),
    (r"pow\((.*?), *-0.30*e1\)",                  r"0.1e1/POW_3(\1)"),
    (r"pow\((.*?), *0.50*e0\)",                   r"sqrt(\1)"),
    (r"pow\((.*?), *-0.50*e0\)",                  r"0.1e1/sqrt(\1)"),
    (r"pow\((.*?), *0.1e1 \/ 0.3e1\)",            r"POW_1_3(\1)"),
    (r"pow\((.*?), *-0.1e1 \/ 0.3e1\)",           r"0.1e1/POW_1_3(\1)"),
    (r"pow\((.*?), *0.1e1 \/ 0.4e1\)",            r"POW_1_4(\1)"),
    (r"pow\((.*?), *-0.1e1 \/ 0.4e1\)",           r"0.1e1/POW_1_4(\1)"),
    (r"pow\((.*?), *0.666666666666666666.e0\)",   r"POW_2_3(\1)"),
    (r"pow\((.*?), *-0.666666666666666666.e0\)",  r"0.1e1 / POW_2_3(\1)"),
    (r"pow\((.*?), *0.333333333333333333.e0\)",   r"POW_1_3($\)"),
    (r"pow\((.*?), *-0.333333333333333333.e0\)",  r"0.1e1 / POW_1_3(\1)"),
    (r"pow\((.*?), *0.1333333333333333333.e1\)",  r"POW_4_3(\1)"),
    (r"pow\((.*?), *-0.1333333333333333333.e1\)", r"0.1e1 / POW_4_3(\1)"),
    (r"pow\((.*?), *0.1666666666666666666.e1\)",  r"POW_5_3(\1)"),
    (r"pow\((.*?), *-0.1666666666666666666.e1\)", r"0.1e1 / POW_5_3(\1)"),
    (r"pow\((.*?), *0.2333333333333333333.e1\)",  r"POW_7_3(\1)"),
    (r"pow\((.*?), *-0.2333333333333333333.e1\)", r"0.1e1 / POW_7_3(\1)"),
    # cleaning up constant expressions
    (r"sqrt\(0.2e1\)",            r"M_SQRT2"),
    (r"POW_1_3\(0.2e1\)",         r"M_CBRT2"),
    (r"POW_1_3\(0.3e1\)",         r"M_CBRT3"),
    (r"POW_1_3\(0.4e1\)",         r"M_CBRT4"),
    (r"POW_1_3\(0.5e1\)",         r"M_CBRT5"),
    (r"POW_1_3\(0.6e1\)",         r"M_CBRT6"),
    (r"POW_1_3\(M_PI\)",          r"M_CBRTPI"),
  )

  # zk_0_ unfortunatly appears in some expressions
  res = re.search(r"zk_0_ = (.*);", text)
  if res:
    text = re.sub(r"zk_0_(?! =)", "(" + res.group(1) + ")", text)

  # standard replacements
  for str1, str2 in math_replace:
    text = re.sub(str1, str2, text)

  # other specific replacements
  for str1, str2 in extra_replace:
    text = re.sub(str1, str2, text)
  
  return text


def maple_run(params, mtype, code, derivatives, start_order):
  '''Creates the maple file, runs maple, and returns the definition
  of the variables and the c-code
  '''
  
  # open maple file
  from tempfile import mkstemp
  fd, mfilename = mkstemp(suffix=".mpl", text=True)
  fh = os.fdopen(fd, "w")
  
  fh.write('''
Polarization := "{}":
Digits := 20:             (* constants will have 20 digits *)
interface(warnlevel=0):   (* supress all warnings          *)
with(CodeGeneration):

$include <{}.mpl>

{}
'''.format(mtype, params["functional"], code))
  fh.close()

  # include dirs for maple
  incdirs = ("maple",
             "maple/lda_exc",  "maple/lda_vxc",
             "maple/gga_exc",  "maple/gga_vxc",
             "maple/mgga_exc", "maple/mgga_vxc"
  )
  maple_inc = ["-I" + params["srcdir"] + "/" + i for i in incdirs]

  # run maple
  run = subprocess.run(
    ["maple"] + maple_inc + ["-q", "-u", mfilename], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
  os.remove(mfilename)
  c_code = run.stdout

  test_1 = ("zk", "vrho", "v2rho2", "v3rho3", "v4rho4", "v5rho5")
  total_order = start_order

  variables  = ["", "", "", "", "", ""]
  n_var = [0, 0, 0, 0, 0, 0]
  new_c_code = ["", "", "", "", "", ""]

  # this adds a new definition of a local variable
  def add_variable(to_add):
    # define 8 variables per line
    if n_var[total_order] % 8 == 0:
      if n_var[total_order] != 0:
        variables[total_order] += ";\n"
      variables[total_order] += "  double "
    else:
      variables[total_order] += ", "
    n_var[total_order] += 1

    variables[total_order] += to_add
  
  # for avoiding compilation when high order derivatives are enabled
  #if start_order != 0:
  #  new_c_code[total_order] += "  if(order < " + str(start_order) + ") return;\n\n\n"

  # we check for strings like 'vrho_0_ = ' and put some
  # relevant if conditions in front
  for line in c_code.splitlines():

    found = False
    # for each order
    for der_order in derivatives:
      # Search the last derivative for each order
      last_derivative = der_order[-1][1]
      # for unpolarized calculation, last derivative is the '0'
      if mtype != "pol":
        last_derivative = re.sub(r"_\d+_", "_0_", last_derivative)

      new_order = re.match(r"\s*" + last_derivative + r"\s*=", line) is not None

      # for each of the derivatives in a given order
      for der in der_order:
        varname  = re.sub(r"_.*", "", der[1])
        varorder = re.sub(r".*_(\d+)_", r"\1", der[1])
        
        # search for a vrho = statement
        if re.match(r"\s*?" + der[1] + r"\s*=", line):

          if (mtype == "pol") or (varorder == "0"):
            # we define a new variable (such as tvrho0) to keep the value
            add_variable("t" + varname + varorder)
            line = re.sub(r"(\S+)_(\d+)_\s*=\s*(.*);",
                          "t" + varname + varorder + r" = \3;", line)
            new_c_code[total_order] += "  " + line + "\n\n"

            # build the if clause to assign the variable
            test = "out->" + varname + " != NULL"

            if not re.search(r"lapl", der[1]) is None:
              test += " && (p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)"

            if not re.search(r"tau", der[1]) is None:
              test += " && (p->info->flags & XC_FLAGS_NEEDS_TAU)"

            test += " && (p->info->flags & XC_FLAGS_HAVE_" + \
              der_name[total_order] + ")"
            new_c_code[total_order] += "  if(" + test + ")\n"

            # add instead of assigning. We are still missing a global constant
            # that can be useful in building hybrid combinations
            new_c_code[total_order] += "    out->{}[ip*p->dim.{} + {}] += t{}{};\n\n".format(varname, varname, varorder, varname, varorder)
            
          found = True
          break

        # find if vrho_0_ is on the right of the = sign
        # this has necessarily to be defined before the
        # left-hand side of the assignement
        while re.search(r"=.*" + varname + r"_\d+_", line):
          line = re.sub(r"(=.*)" + varname + r"_(\d+)_", r"\1t" + varname + r"\2", line)

      # if last variable of this order increment total_order
      if new_order:
        variables[total_order] += ";\n"
        total_order += 1
        
    if not found:
      res = re.match(r"(t\d+) =", line)
      if res: add_variable(res.group(1))

      new_c_code[total_order] += "  " + line + "\n"

  # perform the necessary replacements
  for i in range(len(new_c_code)):
    new_c_code[i] = maple2c_replace(new_c_code[i], params["replace"])

  return variables, new_c_code


def maple2c_run(params, variables, derivatives, variants, start_order, input_args, output_args):

  # open file to write to
  fname = params['srcdir'] + "/src/maple2c/" + \
    params['functype']  + "/" + params['functional'] + ".c"

  out = open(fname, "w")
  if out is None:
    print("Could not open file '" + fname + "' for writing")
    sys.exit(1)
    
  print_c_header(params, out)

  test_2 = ("EXC", "VXC", "FXC", "KXC", "LXC", "MXC")

  out.write("#define MAPLE2C_FLAGS (")
  for i in range(start_order, params['maxorder'] + 1):
    if i != start_order:
      out.write(" | ")
    out.write("XC_FLAGS_I_HAVE_" + test_2[i])
  out.write(")\n\n")

  for mtype, code in variants.items():
    vars_def, c_code = maple_run(params, mtype, code, derivatives, start_order)

    for order in range(start_order, params['maxorder'] + 1):
      out.write('''
#ifndef XC_DONT_COMPILE_{}
GPU_DEVICE_FUNCTION static inline void
func_{}_{}(const xc_func_type *p, size_t ip, {}, {})
{{
'''.format(der_name[order].upper(),
           der_name[order].lower(), mtype,
           input_args, output_args))

      # we first print the declaration of the derivatives
      for order2 in range(start_order, order + 1):
         out.write(vars_def[order2] + "\n")
      
      # we now print the prefix defined in the .mpl code
      out.write(params["prefix"] + "\n")

      # and now the c_code
      for order2 in range(start_order, order + 1):
         out.write(c_code[order2])

      out.write("}\n\n")
      out.write("#endif\n\n")
      

  out.close()
